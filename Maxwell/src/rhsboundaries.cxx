#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Maxwell {
using namespace Loop;

extern "C" void Maxwell_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_RHSBoundaries;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 1, 1> dtdyz_(cctkGH, dtdyz);
  const GF3D<CCTK_REAL, 1, 0, 1> dtdzx_(cctkGH, dtdzx);
  const GF3D<CCTK_REAL, 1, 1, 0> dtdxy_(cctkGH, dtdxy);

  const GF3D<CCTK_REAL, 0, 1, 1> dtbyz_(cctkGH, dtbyz);
  const GF3D<CCTK_REAL, 1, 0, 1> dtbzx_(cctkGH, dtbzx);
  const GF3D<CCTK_REAL, 1, 1, 0> dtbxy_(cctkGH, dtbxy);

  loop_bnd<0, 1, 1>(cctkGH, [&](const PointDesc &p) { dtdyz_(p.I) = 0; });
  loop_bnd<1, 0, 1>(cctkGH, [&](const PointDesc &p) { dtdzx_(p.I) = 0; });
  loop_bnd<1, 1, 0>(cctkGH, [&](const PointDesc &p) { dtdxy_(p.I) = 0; });
  loop_bnd<0, 1, 1>(cctkGH, [&](const PointDesc &p) { dtbyz_(p.I) = 0; });
  loop_bnd<1, 0, 1>(cctkGH, [&](const PointDesc &p) { dtbzx_(p.I) = 0; });
  loop_bnd<1, 1, 0>(cctkGH, [&](const PointDesc &p) { dtbxy_(p.I) = 0; });
}

} // namespace Maxwell
