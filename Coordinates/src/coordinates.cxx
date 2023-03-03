#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Coordinates {

extern "C" void Coordinates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Coordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> vcoordx_(cctkGH, vcoordx);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> vcoordy_(cctkGH, vcoordy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> vcoordz_(cctkGH, vcoordz);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> ccoordx_(cctkGH, ccoordx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> ccoordy_(cctkGH, ccoordy);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> ccoordz_(cctkGH, ccoordz);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> cvol_(cctkGH, cvol);

  Loop::loop_all<0, 0, 0>(
      cctkGH, [&](const Loop::PointDesc &p) { vcoordx_(p.I) = p.x; });
  Loop::loop_all<0, 0, 0>(
      cctkGH, [&](const Loop::PointDesc &p) { vcoordy_(p.I) = p.y; });
  Loop::loop_all<0, 0, 0>(
      cctkGH, [&](const Loop::PointDesc &p) { vcoordz_(p.I) = p.z; });

  Loop::loop_all<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { ccoordx_(p.I) = p.x; });
  Loop::loop_all<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { ccoordy_(p.I) = p.y; });
  Loop::loop_all<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { ccoordz_(p.I) = p.z; });

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    cvol_(p.I) = p.dx * p.dy * p.dz;
  });
}

} // namespace Coordinates
