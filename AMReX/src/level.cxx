#include "schedule.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop.hxx>

namespace AMReX {
using namespace std;

extern "C" void AMReX_InitError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  Loop::loop_all<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { regrid_error[p.idx] = 0; });
}

extern "C" void AMReX_SetLevel(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int levfac = cctk_levfac[0];
  int lev = 0;
  while (levfac > 1) {
    levfac >>= 1;
    lev += 1;
  }

  Loop::loop_all<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { refinement_level[p.idx] = lev; });
}

} // namespace AMReX
