#include <../../CarpetX/src/driver.hxx>
#include <../../CarpetX/src/schedule.hxx>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments_Checked.h>
#include <util_Table.h>

#include <AMReX_MultiFab.H>

#include <cassert>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace ODESolver {
using namespace amrex;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// A state vector component, with mfabs for each level, group, and variable
struct statecomp_t {

  vector<string> groupnames;
  vector<int> groupids;
  vector<MultiFab *> mfabs;

private:
  // These will be automatically freed when this object is deallocated
  vector<unique_ptr<MultiFab> > owned_stuff;

public:
  void check_valid() const;

  statecomp_t copy() const;

  static void axpy(const statecomp_t &y, CCTK_REAL alpha, const statecomp_t &x);
  static void lincomb(const statecomp_t &z, CCTK_REAL alpha,
                      const statecomp_t &x, CCTK_REAL beta,
                      const statecomp_t &y);
};

////////////////////////////////////////////////////////////////////////////////

// Ensure a state vector has valid data everywhere
void statecomp_t::check_valid() const {
  for (const int groupid : groupids) {
    if (groupid >= 0) {
      for (const auto &leveldata : CarpetX::ghext->leveldata) {
        const auto &groupdata = *leveldata.groupdata.at(groupid);
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          const int tl = 0;
          CarpetX::check_valid(leveldata, groupdata, vi, tl, [&]() {
            return "ODESolver before calculating state vector";
          });
        }
      }
    }
  }
}

// Copy state vector into newly allocate memory
statecomp_t statecomp_t::copy() const {
  check_valid();
  const size_t size = mfabs.size();
  statecomp_t result;
  result.groupnames = groupnames;
  result.groupids.resize(size, -1);
  result.mfabs.reserve(size);
  result.owned_stuff.reserve(size);
  for (size_t n = 0; n < size; ++n) {
    const auto &x = mfabs.at(n);
    if (x->contains_nan())
      CCTK_VERROR("statecomp_t::copy.x: Group %s contains nans",
                  groupnames.at(n).c_str());
    auto y = make_unique<MultiFab>(x->boxArray(), x->DistributionMap(),
                                   x->nComp(), x->nGrowVect());
    MultiFab::Copy(*y, *x, 0, 0, y->nComp(), y->nGrowVect());
    if (y->contains_nan())
      CCTK_VERROR("statecomp_t::copy.y: Group %s contains nans",
                  result.groupnames.at(n).c_str());
    result.mfabs.push_back(y.get());
    result.owned_stuff.push_back(move(y));
  }
  result.check_valid();
  return result;
}

// y += alpha * x
void statecomp_t::axpy(const statecomp_t &y, const CCTK_REAL alpha,
                       const statecomp_t &x) {
  x.check_valid();
  y.check_valid();
  const size_t size = y.mfabs.size();
  assert(x.mfabs.size() == size);
  for (size_t n = 0; n < size; ++n) {
    assert(x.mfabs.at(n)->nGrowVect() == y.mfabs.at(n)->nGrowVect());
    if (x.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::axpy.x: Group %s contains nans",
                  x.groupnames.at(n).c_str());
    MultiFab::Saxpy(*y.mfabs.at(n), alpha, *x.mfabs.at(n), 0, 0,
                    y.mfabs.at(n)->nComp(), y.mfabs.at(n)->nGrowVect());
    if (y.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::axpy.y: Group %s contains nans",
                  y.groupnames.at(n).c_str());
  }
  y.check_valid();
}

// z = alpha * x + beta * y
void statecomp_t::lincomb(const statecomp_t &z, const CCTK_REAL alpha,
                          const statecomp_t &x, const CCTK_REAL beta,
                          const statecomp_t &y) {
  x.check_valid();
  y.check_valid();
  const size_t size = z.mfabs.size();
  assert(x.mfabs.size() == size);
  assert(y.mfabs.size() == size);
  for (size_t n = 0; n < size; ++n) {
    assert(x.mfabs.at(n)->nGrowVect() == z.mfabs.at(n)->nGrowVect());
    assert(y.mfabs.at(n)->nGrowVect() == z.mfabs.at(n)->nGrowVect());
    if (x.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::lincomb.x: Group %s contains nans",
                  x.groupnames.at(n).c_str());
    if (y.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::lincomb.y: Group %s contains nans",
                  y.groupnames.at(n).c_str());
    MultiFab::LinComb(*z.mfabs.at(n), alpha, *x.mfabs.at(n), 0, beta,
                      *y.mfabs.at(n), 0, 0, z.mfabs.at(n)->nComp(),
                      z.mfabs.at(n)->nGrowVect());
    if (z.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::lincomb.z: Group %s contains nans",
                  z.groupnames.at(n).c_str());
  }
  z.check_valid();
}

////////////////////////////////////////////////////////////////////////////////

int get_group_rhs(const int gi) {
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  vector<char> rhs_buf(1000);
  const int iret =
      Util_TableGetString(tags, rhs_buf.size(), rhs_buf.data(), "rhs");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    rhs_buf[0] = '\0'; // default: empty (no RHS)
  } else if (iret >= 0) {
    // do nothing
  } else {
    assert(0);
  }

  const string str(rhs_buf.data());
  if (str.empty())
    return -1; // No RHS specified

  auto str1 = str;
  if (str1.find(':') == string::npos) {
    const char *impl = CCTK_GroupImplementationI(gi);
    str1 = string(impl) + "::" + str1;
  }
  const int gi1 = CCTK_GroupIndex(str1.c_str());
  assert(gi1 >= 0); // Check fluxes are valid groups
  const int flux = gi1;

  assert(flux != gi);

  return flux;
}

///////////////////////////////////////////////////////////////////////////////

extern "C" void ODESolvers_Solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ODESolvers_Solve;
  DECLARE_CCTK_PARAMETERS;

  static bool did_output = false;
  if (verbose || !did_output)
    CCTK_VINFO("Integrator is %s", method);
  did_output = true;

  const CCTK_REAL dt = cctk_delta_time;
  const int tl = 0;

  statecomp_t var, rhs;
  int nvars = 0;
  for (const auto &leveldata : CarpetX::ghext->leveldata) {
    for (const auto &groupdataptr : leveldata.groupdata) {
      // TODO: add support for evolving grid scalars
      if (groupdataptr == nullptr)
        continue;

      const auto &groupdata = *groupdataptr;
      const int rhs_gi = get_group_rhs(groupdata.groupindex);
      if (rhs_gi >= 0) {
        const auto &rhs_groupdata = *leveldata.groupdata.at(rhs_gi);
        assert(rhs_groupdata.numvars == groupdata.numvars);
        var.groupnames.push_back(CCTK_GroupName(groupdata.groupindex));
        var.groupids.push_back(groupdata.groupindex);
        var.mfabs.push_back(groupdata.mfab.at(tl).get());
        rhs.groupnames.push_back(CCTK_GroupName(rhs_groupdata.groupindex));
        rhs.groupids.push_back(rhs_groupdata.groupindex);
        rhs.mfabs.push_back(rhs_groupdata.mfab.at(tl).get());
        if (leveldata.level == 0)
          nvars += groupdata.numvars;
      }
    }
  }
  if (verbose)
    CCTK_VINFO("  Integrating %d variables", nvars);
  if (nvars == 0)
    CCTK_VWARN(CCTK_WARN_ALERT, "Integrating %d variables", nvars);

  const CCTK_REAL saved_time = cctkGH->cctk_time;
  const CCTK_REAL old_time = cctkGH->cctk_time - dt;
  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time;
  // Calculate first RHS
  if (verbose)
    CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
  CallScheduleGroup(cctkGH, "ODESolvers_RHS");

  if (CCTK_EQUALS(method, "constant")) {

    // do nothing

  } else if (CCTK_EQUALS(method, "Euler")) {

    // Add scaled RHS to state vector
    statecomp_t::axpy(var, dt, rhs);

  } else if (CCTK_EQUALS(method, "RK2")) {

    const auto old = var.copy();

    // Add scaled RHS to state vector
    statecomp_t::axpy(var, dt / 2, rhs);
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) += dt / 2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Calculate new state vector
    statecomp_t::lincomb(var, 1, old, dt, rhs);

  } else if (CCTK_EQUALS(method, "RK4")) {

    const auto old = var.copy();
    const auto k1 = rhs.copy();

    // Step 2

    // Add scaled RHS to state vector
    statecomp_t::axpy(var, dt / 2, rhs);
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto k2 = rhs.copy();

    // Step 3

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 1, old, dt / 2, rhs);
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt / 2;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #3 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto k3 = rhs.copy();

    // Step 4

    // Add scaled RHS to state vector
    statecomp_t::lincomb(var, 1, old, dt, rhs);
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = old_time + dt;
    CallScheduleGroup(cctkGH, "ODESolvers_PostStep");

    if (verbose)
      CCTK_VINFO("Calculating RHS #4 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    const auto &k4 = rhs;

    // Calculate new state vector
    statecomp_t::lincomb(var, 1, old, dt / 6, k1);
    statecomp_t::axpy(var, dt / 3, k2);
    statecomp_t::axpy(var, dt / 3, k3);
    statecomp_t::axpy(var, dt / 6, k4);

  } else {
    assert(0);
  }

  // Reset current time
  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = saved_time;
  // Apply last boundary conditions
  CallScheduleGroup(cctkGH, "ODESolvers_PostStep");
  if (verbose)
    CCTK_VINFO("Calculated new state at t=%g", double(cctkGH->cctk_time));

  // TODO: Update time here, and not during time level cycling in the driver
}

} // namespace ODESolver
