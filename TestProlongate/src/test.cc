#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>

namespace TestProlongate {
using namespace std;

extern "C" void TestProlongate_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  // get prolongation order from driver, the parmeter is private since really
  // there is normally no reason to depend on it
  int order_type;
  const void *order_p = CCTK_ParameterGet("prolongation_order", "AMReX",
                                          &order_type);
  assert(order_p);
  assert(order_type == PARAMETER_INT);
  const CCTK_INT operator_order = *static_cast<const CCTK_INT*>(order_p);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    data[p.idx] = pow(p.x * p.y * p.z, operator_order);
  });

  *max_diff = 0.;
}

extern "C" void TestProlongate_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  return; // do nothing
}
extern "C" void TestProlongate_Regrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (cctk_iteration < regrid_after)
    return;

  CCTK_VINFO("Setting grid at %d\n", cctk_iteration);
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if (fabs(p.x) <= refined_radius && fabs(p.y) <= refined_radius &&
        fabs(p.z) <= refined_radius) {
      regrid_error[p.idx] = 1e3;
    } else {
      regrid_error[p.idx] = 0.;
    }
  });
}

extern "C" void TestProlongate_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  // there is no reduction yet, so for now I require that only 1 MPI rank be in
  // use
  assert(CCTK_nProcs(cctkGH) == 1 && "Must be single rank");

  // get prolongation order from driver, the parmeter is private since really
  // there is normally no reason to depend on it
  int order_type;
  const void *order_p = CCTK_ParameterGet("prolongation_order", "AMReX",
                                          &order_type);
  assert(order_p);
  assert(order_type == PARAMETER_INT);
  const CCTK_INT operator_order = *static_cast<const CCTK_INT*>(order_p);

  CCTK_REAL my_max_diff = 0;
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const CCTK_REAL diff = fabs(data[p.idx] -
                                pow(p.x * p.y * p.z, operator_order));
    my_max_diff = max(diff, my_max_diff);
  });

  #pragma omp critical
  {
    *max_diff = max(*max_diff, my_max_diff);
  }
}

} // namespace TestProlongate
