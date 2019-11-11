#include <cctk.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>

namespace TestODESolvers {
using namespace std;

extern "C" void TestODESolvers_initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers_initial;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL u0 = pow(1 + cctk_time, order);

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
        state[ind] = u0;
      }
    }
  }
}

extern "C" void TestODESolvers_boundary(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers_boundary;
  DECLARE_CCTK_PARAMETERS;

  // do nothing
}

extern "C" void TestODESolvers_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers_RHS;
  DECLARE_CCTK_PARAMETERS;

  // u(t) = (1+t)^p
  // d/dt u = p (1+t)^(p-1)

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
        rhs[ind] = order * pow(1 + cctk_time, order - 1);
      }
    }
  }
}

extern "C" void TestODESolvers_error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TestODESolvers_error;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL u0 = pow(1 + cctk_time, order);

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        int ind = CCTK_GFINDEX3D(cctkGH, i, j, k);
        error[ind] = state[ind] - u0;
      }
    }
  }
}

} // namespace TestODESolvers
