#include "schedule.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace AMReX {
using namespace std;

extern "C" void AMReX_SetLevel(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int levfac = cctk_levfac[0];
  int lev = 0;
  while (levfac > 1) {
    levfac >>= 1;
    lev += 1;
  }

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        refinement_level[idx] = lev;
      }
    }
  }
}

} // namespace AMReX
