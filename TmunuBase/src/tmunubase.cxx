#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace TmunuBase {
using namespace Loop;

extern "C" void TmunuBase_ZeroTmunu(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_TmunuBase_ZeroTmunu;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> eTtt_(cctkGH, eTtt);
  const GF3D<CCTK_REAL, 0, 0, 0> eTtx_(cctkGH, eTtx);
  const GF3D<CCTK_REAL, 0, 0, 0> eTty_(cctkGH, eTty);
  const GF3D<CCTK_REAL, 0, 0, 0> eTtz_(cctkGH, eTtz);
  const GF3D<CCTK_REAL, 0, 0, 0> eTxx_(cctkGH, eTxx);
  const GF3D<CCTK_REAL, 0, 0, 0> eTxy_(cctkGH, eTxy);
  const GF3D<CCTK_REAL, 0, 0, 0> eTxz_(cctkGH, eTxz);
  const GF3D<CCTK_REAL, 0, 0, 0> eTyy_(cctkGH, eTyy);
  const GF3D<CCTK_REAL, 0, 0, 0> eTyz_(cctkGH, eTyz);
  const GF3D<CCTK_REAL, 0, 0, 0> eTzz_(cctkGH, eTzz);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTtt_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTtx_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTty_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTtz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTxx_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTxy_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTxz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTyy_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTyz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { eTzz_(p.I) = 0; });
}

} // namespace TmunuBase
