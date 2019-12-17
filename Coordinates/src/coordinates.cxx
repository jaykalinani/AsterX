#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace Coordinates {

extern "C" void Coordinates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Coordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> coordx_(cctkGH, coordx);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> coordy_(cctkGH, coordy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> coordz_(cctkGH, coordz);

  Loop::loop_all<0, 0, 0>(
      cctkGH, [&](const Loop::PointDesc &p) { coordx_(p.I) = p.x; });
  Loop::loop_all<0, 0, 0>(
      cctkGH, [&](const Loop::PointDesc &p) { coordy_(p.I) = p.y; });
  Loop::loop_all<0, 0, 0>(
      cctkGH, [&](const Loop::PointDesc &p) { coordz_(p.I) = p.z; });
}

} // namespace Coordinates
