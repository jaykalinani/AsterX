#include <../../CarpetX/src/driver.hxx>

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
  vector<string> varnames;
  vector<MultiFab *> mfabs;

private:
  // These will be automatically freed when this object is deallocated
  vector<unique_ptr<MultiFab> > owned_stuff;

public:
  statecomp_t copy() const;

  static void axpy(const statecomp_t &y, CCTK_REAL alpha, const statecomp_t &x);
  static void lincomb(const statecomp_t &z, CCTK_REAL alpha,
                      const statecomp_t &x, CCTK_REAL beta,
                      const statecomp_t &y);
};

////////////////////////////////////////////////////////////////////////////////

// Copy state vector into newly allocate memory
statecomp_t statecomp_t::copy() const {
  const size_t size = mfabs.size();
  statecomp_t result;
  result.varnames = varnames;
  result.mfabs.reserve(size);
  result.owned_stuff.reserve(size);
  for (size_t n = 0; n < size; ++n) {
    const auto &x = mfabs.at(n);
    if (x->contains_nan())
      CCTK_VERROR("statecomp_t::copy.x: Variable %s contains nans",
                  varnames.at(n).c_str());
    auto y = make_unique<MultiFab>(x->boxArray(), x->DistributionMap(),
                                   x->nComp(), x->nGrowVect());
    MultiFab::Copy(*y, *x, 0, 0, y->nComp(), y->nGrowVect());
    if (y->contains_nan())
      CCTK_VERROR("statecomp_t::copy.y: Variable %s contains nans",
                  result.varnames.at(n).c_str());
    result.mfabs.push_back(y.get());
    result.owned_stuff.push_back(move(y));
  }
  return result;
}

// y += alpha * x
void statecomp_t::axpy(const statecomp_t &y, CCTK_REAL alpha,
                       const statecomp_t &x) {
  const size_t size = y.mfabs.size();
  assert(x.mfabs.size() == size);
  for (size_t n = 0; n < size; ++n) {
    assert(x.mfabs.at(n)->nGrowVect() == y.mfabs.at(n)->nGrowVect());
    if (x.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::axpy.x: Variable %s contains nans",
                  x.varnames.at(n).c_str());
    MultiFab::Saxpy(*y.mfabs.at(n), alpha, *x.mfabs.at(n), 0, 0,
                    y.mfabs.at(n)->nComp(), y.mfabs.at(n)->nGrowVect());
    if (y.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::axpy.y: Variable %s contains nans",
                  y.varnames.at(n).c_str());
  }
}

// z = alpha * x + beta * y
void statecomp_t::lincomb(const statecomp_t &z, CCTK_REAL alpha,
                          const statecomp_t &x, CCTK_REAL beta,
                          const statecomp_t &y) {
  const size_t size = z.mfabs.size();
  assert(x.mfabs.size() == size);
  assert(y.mfabs.size() == size);
  for (size_t n = 0; n < size; ++n) {
    assert(x.mfabs.at(n)->nGrowVect() == z.mfabs.at(n)->nGrowVect());
    assert(y.mfabs.at(n)->nGrowVect() == z.mfabs.at(n)->nGrowVect());
    if (x.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::lincomb.x: Variable %s contains nans",
                  x.varnames.at(n).c_str());
    if (y.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::lincomb.y: Variable %s contains nans",
                  y.varnames.at(n).c_str());
    MultiFab::LinComb(*z.mfabs.at(n), alpha, *x.mfabs.at(n), 0, beta,
                      *y.mfabs.at(n), 0, 0, z.mfabs.at(n)->nComp(),
                      z.mfabs.at(n)->nGrowVect());
    if (z.mfabs.at(n)->contains_nan())
      CCTK_VERROR("statecomp_t::lincomb.z %s contains nans",
                  z.varnames.at(n).c_str());
  }
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

  CCTK_VINFO("Integrator is %s", method);

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const int tl = 0;

  statecomp_t var, rhs;
  for (const auto &leveldata : CarpetX::ghext->leveldata) {
    for (const auto &groupdata : leveldata.groupdata) {
      const int rhs_gi = get_group_rhs(groupdata->groupindex);
      if (rhs_gi >= 0) {
        const auto &rhs_groupdata = leveldata.groupdata.at(rhs_gi);
        assert(rhs_groupdata->numvars == groupdata->numvars);
        for (int vi = 0; vi < groupdata->numvars; ++vi) {
          var.varnames.push_back(
              CCTK_FullVarName(groupdata->firstvarindex + vi));
          var.mfabs.push_back(groupdata->mfab.at(tl).get());
          rhs.varnames.push_back(
              CCTK_FullVarName(rhs_groupdata->firstvarindex + vi));
          rhs.mfabs.push_back(rhs_groupdata->mfab.at(tl).get());
        }
      }
    }
  }
  const int nvars = int(var.mfabs.size());
  CCTK_VINFO("  Integrating %d variables", nvars);
  if (nvars == 0)
    CCTK_VWARN(CCTK_WARN_ALERT, "Integrating %d variables", nvars);

  const CCTK_REAL saved_time = cctkGH->cctk_time;
  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) -= dt;

  if (CCTK_EQUALS(method, "constant")) {

    // do nothing

  } else if (CCTK_EQUALS(method, "Euler")) {

    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Add scaled RHS to state vector
    statecomp_t::axpy(var, dt, rhs);

  } else if (CCTK_EQUALS(method, "RK2")) {

    const auto old = var.copy();

    CCTK_VINFO("Calculating RHS #1 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Add scaled RHS to state vector
    statecomp_t::axpy(var, dt / 2, rhs);
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) += dt / 2;

    CCTK_VINFO("Calculating RHS #2 at t=%g", double(cctkGH->cctk_time));
    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Calculate new state vector
    statecomp_t::lincomb(var, 1, old, dt, rhs);

  } else {
    assert(0);
  }

  // Reset current time
  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = saved_time;
  CCTK_VINFO("Calculated new state at t=%g", double(cctkGH->cctk_time));

  // TODO: Update time here, and not during time level cycling in the driver
}

} // namespace ODESolver
