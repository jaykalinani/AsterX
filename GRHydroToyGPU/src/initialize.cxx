#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

extern "C" void GRHydroToyGPU_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout(cctkGH, cell_centred);

  const GF3D2<CCTK_REAL> gf_dens(gf_layout, dens);
  const GF3D2<CCTK_REAL> gf_momx(gf_layout, momx);
  const GF3D2<CCTK_REAL> gf_momy(gf_layout, momy);
  const GF3D2<CCTK_REAL> gf_momz(gf_layout, momz);
  const GF3D2<CCTK_REAL> gf_tau(gf_layout, tau);

  if (CCTK_EQUALS(setup, "equilibrium")) {

    grid.loop_int_device<1, 1, 1>(grid.nghostzones,
                                  [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
                                      CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                        gf_dens(p.I) = 1.0;
                                        gf_momx(p.I) = 0.0;
                                        gf_momy(p.I) = 0.0;
                                        gf_momz(p.I) = 0.0;
                                        gf_tau(p.I) = 1.0;
                                      });

  } else if (CCTK_EQUALS(setup, "sound wave")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                              const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gf_dens(p.I) = 1.0;
          gf_momx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          gf_momy(p.I) = 0.0;
          gf_momz(p.I) = 0.0;
          gf_tau(p.I) = 1.0; // should add kinetic energy here
        });

  } else if (CCTK_EQUALS(setup, "shock tube")) {

    grid.loop_int_device<1, 1, 1>(grid.nghostzones,
                                  [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
                                      CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                        if (p.x <= 0.0) {
                                          gf_dens(p.I) = 2.0;
                                          gf_momx(p.I) = 0.0;
                                          gf_momy(p.I) = 0.0;
                                          gf_momz(p.I) = 0.0;
                                          gf_tau(p.I) = 2.0;
                                        } else {
                                          gf_dens(p.I) = 1.0;
                                          gf_momx(p.I) = 0.0;
                                          gf_momy(p.I) = 0.0;
                                          gf_momz(p.I) = 0.0;
                                          gf_tau(p.I) = 1.0;
                                        }
                                      });

  } else if (CCTK_EQUALS(setup, "spherical shock")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                              const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL r2 = pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2);
          if (r2 <= pow(shock_radius, 2)) {
            gf_dens(p.I) = 2.0;
            gf_momx(p.I) = 0.0;
            gf_momy(p.I) = 0.0;
            gf_momz(p.I) = 0.0;
            gf_tau(p.I) = 2.0;
          } else {
            gf_dens(p.I) = 1.0;
            gf_momx(p.I) = 0.0;
            gf_momy(p.I) = 0.0;
            gf_momz(p.I) = 0.0;
            gf_tau(p.I) = 1.0;
          }
        });

  } else {
    CCTK_ERROR("Internal error");
  }
}

} // namespace GRHydroToyGPU
