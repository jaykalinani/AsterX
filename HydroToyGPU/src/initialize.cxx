#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace HydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

extern "C" void HydroToyGPU_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout(cctkGH, cell_centred);

  const GF3D2<CCTK_REAL> gf_rho(gf_layout, rho);
  const GF3D2<CCTK_REAL> gf_momx(gf_layout, momx);
  const GF3D2<CCTK_REAL> gf_momy(gf_layout, momy);
  const GF3D2<CCTK_REAL> gf_momz(gf_layout, momz);
  const GF3D2<CCTK_REAL> gf_etot(gf_layout, etot);

  if (CCTK_EQUALS(setup, "equilibrium")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=](const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE {
              gf_rho(p.I) = 1.0;
              gf_momx(p.I) = 0.0;
              gf_momy(p.I) = 0.0;
              gf_momz(p.I) = 0.0;
              gf_etot(p.I) = 1.0;
            });

  } else if (CCTK_EQUALS(setup, "sound wave")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=](const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE {
              gf_rho(p.I) = 1.0;
              gf_momx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
              gf_momy(p.I) = 0.0;
              gf_momz(p.I) = 0.0;
              gf_etot(p.I) = 1.0; // should add kinetic energy here
            });

  } else if (CCTK_EQUALS(setup, "shock tube")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=](const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE {
              if (p.x <= 0.0) {
                gf_rho(p.I) = 2.0;
                gf_momx(p.I) = 0.0;
                gf_momy(p.I) = 0.0;
                gf_momz(p.I) = 0.0;
                gf_etot(p.I) = 2.0;
              } else {
                gf_rho(p.I) = 1.0;
                gf_momx(p.I) = 0.0;
                gf_momy(p.I) = 0.0;
                gf_momz(p.I) = 0.0;
                gf_etot(p.I) = 1.0;
              }
            });

  } else if (CCTK_EQUALS(setup, "spherical shock")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=](const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE {
              CCTK_REAL r2 = pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2);
              if (r2 <= pow(shock_radius, 2)) {
                gf_rho(p.I) = 2.0;
                gf_momx(p.I) = 0.0;
                gf_momy(p.I) = 0.0;
                gf_momz(p.I) = 0.0;
                gf_etot(p.I) = 2.0;
              } else {
                gf_rho(p.I) = 1.0;
                gf_momx(p.I) = 0.0;
                gf_momy(p.I) = 0.0;
                gf_momz(p.I) = 0.0;
                gf_etot(p.I) = 1.0;
              }
            });

  } else {
    CCTK_ERROR("Internal error");
  }
}

} // namespace HydroToyGPU
