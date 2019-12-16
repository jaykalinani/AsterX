#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace ADMBase {

extern "C" void ADMBase_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_data;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> gxx_(cctkGH, gxx);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> gxy_(cctkGH, gxy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> gxz_(cctkGH, gxz);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> gyy_(cctkGH, gyy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> gyz_(cctkGH, gyz);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> gzz_(cctkGH, gzz);

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> kxx_(cctkGH, kxx);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> kxy_(cctkGH, kxy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> kxz_(cctkGH, kxz);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> kyy_(cctkGH, kyy);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> kyz_(cctkGH, kyz);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> kzz_(cctkGH, kzz);

  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { gxx_(p.I) = 1; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { gxy_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { gxz_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { gyy_(p.I) = 1; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { gyz_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { gzz_(p.I) = 1; });

  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { kxx_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { kxy_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { kxz_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { kyy_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { kyz_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { kzz_(p.I) = 0; });
}

extern "C" void ADMBase_initial_lapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_lapse;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> alp_(cctkGH, alp);

  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { alp_(p.I) = 1; });
}

extern "C" void ADMBase_initial_dtlapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_dtlapse;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> dtalp_(cctkGH, dtalp);

  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { dtalp_(p.I) = 0; });
}

extern "C" void ADMBase_initial_shift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_shift;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> betax_(cctkGH, betax);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> betay_(cctkGH, betay);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> betaz_(cctkGH, betaz);

  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { betax_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { betay_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { betaz_(p.I) = 0; });
}

extern "C" void ADMBase_initial_dtshift(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_ADMBase_initial_dtshift;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> dtbetax_(cctkGH, dtbetax);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> dtbetay_(cctkGH, dtbetay);
  const Loop::GF3D<CCTK_REAL, 0, 0, 0> dtbetaz_(cctkGH, dtbetaz);

  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { dtbetax_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { dtbetay_(p.I) = 0; });
  Loop::loop_all<0, 0, 0>(cctkGH,
                          [&](const Loop::PointDesc &p) { dtbetaz_(p.I) = 0; });
}

} // namespace ADMBase
