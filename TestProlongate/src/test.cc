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

  int prolongation_type_type;
  const void *prolongation_type_p =
    CCTK_ParameterGet("prolongation_type", "AMReX", &prolongation_type_type);
  assert(prolongation_type_p);
  assert(prolongation_type_type == PARAMETER_KEYWORD);
  const char* prolongation_type =
    *static_cast<const char* const*>(prolongation_type_p);
  const bool conservative_prolongation =
    CCTK_EQUALS(prolongation_type, "conservative");

  // the grid stores the average values of the underlying function
  // (x*y*z)**n which amounts to storing differences of the anti-derivative
  // 1/(n+1)**3 * (x*y*z)**(n+1)
  auto avg_fun = [&](CCTK_REAL xi, int dir) {
    return 1./(operator_order + 1) *
           (pow(xi+CCTK_DELTA_SPACE(dir)/2., operator_order + 1) -
            pow(xi-CCTK_DELTA_SPACE(dir)/2., operator_order + 1)) /
           CCTK_DELTA_SPACE(dir);
  };

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    CCTK_REAL good_data;
    if (conservative_prolongation) {
      good_data = avg_fun(p.x, 0) * avg_fun(p.y, 1) * avg_fun(p.z, 2);
    } else {
      good_data = pow(p.x * p.y * p.z, operator_order);
    }
    data[p.idx] = good_data;
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

  CCTK_VINFO("Setting grid at %d", cctk_iteration);
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

  int prolongation_type_type;
  const void *prolongation_type_p =
    CCTK_ParameterGet("prolongation_type", "AMReX", &prolongation_type_type);
  assert(prolongation_type_p);
  assert(prolongation_type_type == PARAMETER_KEYWORD);
  const char* prolongation_type =
    *static_cast<const char* const*>(prolongation_type_p);
  const bool conservative_prolongation =
    CCTK_EQUALS(prolongation_type, "conservative");

  // the grid stores the average values of the underlying function
  // (x*y*z)**n which amounts to storing differences of the anti-derivative
  // 1/(n+1)**3 * (x*y*z)**(n+1)
  auto avg_fun = [&](CCTK_REAL xi, int dir) {
    return 1./(operator_order + 1) *
           (pow(xi+CCTK_DELTA_SPACE(dir)/2., operator_order + 1) -
            pow(xi-CCTK_DELTA_SPACE(dir)/2., operator_order + 1)) /
           CCTK_DELTA_SPACE(dir);
  };

  CCTK_REAL my_max_diff = 0;
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    CCTK_REAL good_data;
    if (conservative_prolongation) {
      good_data = avg_fun(p.x, 0) * avg_fun(p.y, 1) * avg_fun(p.z, 2);
    } else {
      good_data = pow(p.x * p.y * p.z, operator_order);
    }
    const CCTK_REAL diff = fabs(data[p.idx] - good_data);
    my_max_diff = max(diff, my_max_diff);
  });

  #pragma omp critical
  {
    *max_diff = max(*max_diff, my_max_diff);
  }
}

} // namespace TestProlongate
