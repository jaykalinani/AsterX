#include <../../CarpetX/src/driver.hxx>

#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments_Checked.h>

#include <AMReX_MultiFab.H>

#include <cassert>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

namespace ODESolver {
using namespace amrex;
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// A state vector component, with mfabs for each level, group, and variable
struct statecomp_t {
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
  result.mfabs.reserve(size);
  result.owned_stuff.reserve(size);
  for (size_t n = 0; n < size; ++n) {
    const auto &x = mfabs.at(n);
    auto y = make_unique<MultiFab>(x->boxArray(), x->DistributionMap(),
                                   x->nComp(), x->nGrowVect());
    MultiFab::Copy(*y, *x, 0, 0, y->nComp(), y->nGrowVect());
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
    // TODO: Update Saxpy call when using AMReX 19.12 or later
    const auto nghosts = y.mfabs.at(n)->nGrowVect();
    for (int d = 0; d < 3; ++d)
      assert(nghosts[d] == nghosts[0]);
    MultiFab::Saxpy(*y.mfabs.at(n), alpha, *x.mfabs.at(n), 0, 0,
                    y.mfabs.at(n)->nComp(), nghosts[0]);
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
    // TODO: Update Lincomb call when using AMReX 19.12 or later
    const auto nghosts = z.mfabs.at(n)->nGrowVect();
    for (int d = 0; d < 3; ++d)
      assert(nghosts[d] == nghosts[0]);
    MultiFab::LinComb(*z.mfabs.at(n), alpha, *x.mfabs.at(n), 0, beta,
                      *y.mfabs.at(n), 0, 0, z.mfabs.at(n)->nComp(), nghosts[0]);
  }
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void ODESolvers_solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ODESolvers_solve;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const int tl = 0;

  statecomp_t var, rhs;
  for (const auto &leveldata : CarpetX::ghext->leveldata) {
    for (const auto &groupdata : leveldata.groupdata) {
      if (groupdata.groupindex == CCTK_GroupIndex("TestODESolvers::state")) {
        const int rhs_gi = CCTK_GroupIndex("TestODESolvers::rhs");
        const auto &rhs_groupdata = leveldata.groupdata.at(rhs_gi);
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          var.mfabs.push_back(groupdata.mfab.at(tl).get());
          rhs.mfabs.push_back(rhs_groupdata.mfab.at(tl).get());
        }
      }
    }
  }

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

    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Add scaled RHS to state vector
    statecomp_t::axpy(var, dt / 2, rhs);
    *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) += dt / 2;

    CallScheduleGroup(cctkGH, "ODESolvers_RHS");

    // Calculate new state vector
    statecomp_t::lincomb(var, 1, old, dt, rhs);

  } else {
    assert(0);
  }

  // Reset current time
  *const_cast<CCTK_REAL *>(&cctkGH->cctk_time) = saved_time;

  // TODO: Update time here, and not during time level cycling in the driver
}

} // namespace ODESolver
