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

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> regrid_error_(cctkGH, regrid_error);

  Loop::loop<1, 1, 1>(
      cctkGH, where_t::everywhere,
      [&](const Loop::PointDesc &p) { regrid_error_(p.I) = 0; });
}

extern "C" void AMReX_SetLevel(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> refinement_level_(cctkGH,
                                                         refinement_level);

  int levfac = cctk_levfac[0];
  int lev = 0;
  while (levfac > 1) {
    levfac >>= 1;
    lev += 1;
  }

  Loop::loop<1, 1, 1>(
      cctkGH, where_t::everywhere,
      [&](const Loop::PointDesc &p) { refinement_level_(p.I) = lev; });
}

} // namespace AMReX
