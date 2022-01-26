#include <fixmath.hxx>
#include <loop.hxx>

#include "schedule.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace CarpetX {
using namespace std;

extern "C" void CarpetX_InitError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CarpetX_InitError;
  DECLARE_CCTK_PARAMETERS;

  Loop::loop<1, 1, 1>(cctkGH, where_t::everywhere,
                      [&](const Loop::PointDesc &p) { regrid_error(p.I) = 0; });
}

extern "C" void CarpetX_SetLevel(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_CarpetX_SetLevel;
  DECLARE_CCTK_PARAMETERS;

  int levfac = cctk_levfac[0];
  int lev = 0;
  while (levfac > 1) {
    levfac >>= 1;
    lev += 1;
  }

  Loop::loop<1, 1, 1>(
      cctkGH, where_t::everywhere,
      [&](const Loop::PointDesc &p) { refinement_level(p.I) = lev; });
}

} // namespace CarpetX
