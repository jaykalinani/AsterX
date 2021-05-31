#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace ADMBase {
using namespace Loop;

extern "C" void ADMBase_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_data;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> gxx_(cctkGH, gxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gxy_(cctkGH, gxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gxz_(cctkGH, gxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gyy_(cctkGH, gyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gyz_(cctkGH, gyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gzz_(cctkGH, gzz);

  const GF3D<CCTK_REAL, 0, 0, 0> kxx_(cctkGH, kxx);
  const GF3D<CCTK_REAL, 0, 0, 0> kxy_(cctkGH, kxy);
  const GF3D<CCTK_REAL, 0, 0, 0> kxz_(cctkGH, kxz);
  const GF3D<CCTK_REAL, 0, 0, 0> kyy_(cctkGH, kyy);
  const GF3D<CCTK_REAL, 0, 0, 0> kyz_(cctkGH, kyz);
  const GF3D<CCTK_REAL, 0, 0, 0> kzz_(cctkGH, kzz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxx_(p.I) = 1; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxy_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxz_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gyy_(p.I) = 1; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gyz_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gzz_(p.I) = 1; });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxx_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxy_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxz_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kyy_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kyz_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kzz_(p.I) = 0; });
}

extern "C" void ADMBase_initial_lapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_lapse;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> alp_(cctkGH, alp);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { alp_(p.I) = 1; });
}

extern "C" void ADMBase_initial_dtlapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_dtlapse;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> dtalp_(cctkGH, dtalp);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { dtalp_(p.I) = 0; });
}

extern "C" void ADMBase_initial_shift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_shift;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> betax_(cctkGH, betax);
  const GF3D<CCTK_REAL, 0, 0, 0> betay_(cctkGH, betay);
  const GF3D<CCTK_REAL, 0, 0, 0> betaz_(cctkGH, betaz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { betax_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { betay_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { betaz_(p.I) = 0; });
}

extern "C" void ADMBase_initial_dtshift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_dtshift;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> dtbetax_(cctkGH, dtbetax);
  const GF3D<CCTK_REAL, 0, 0, 0> dtbetay_(cctkGH, dtbetay);
  const GF3D<CCTK_REAL, 0, 0, 0> dtbetaz_(cctkGH, dtbetaz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { dtbetax_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { dtbetay_(p.I) = 0; });
  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) { dtbetaz_(p.I) = 0; });
}

} // namespace ADMBase
