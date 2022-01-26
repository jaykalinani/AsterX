#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace HydroInitial {
using namespace std;
using namespace Loop;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void HydroInitial_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroInitial_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout_cell(cctkGH, cell_centred);

  const GF3D2<CCTK_REAL> gf_rho(gf_layout_cell, rho);
  const GF3D2<CCTK_REAL> gf_velx(gf_layout_cell, velx);
  const GF3D2<CCTK_REAL> gf_vely(gf_layout_cell, vely);
  const GF3D2<CCTK_REAL> gf_velz(gf_layout_cell, velz);
  const GF3D2<CCTK_REAL> gf_press(gf_layout_cell, press);
  const GF3D2<CCTK_REAL> gf_eps(gf_layout_cell, eps);

  if (CCTK_EQUALS(initial_hydro, "equilibrium")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gf_rho(p.I) = 1.0;
          gf_velx(p.I) = 0.0;
          gf_vely(p.I) = 0.0;
          gf_velz(p.I) = 0.0;
          gf_press(p.I) = 1.0;
          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I) / (gf_rho(p.I) * (gamma - 1));
        });

  } else if (CCTK_EQUALS(initial_hydro, "sound wave")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gf_rho(p.I) = 1.0;
          gf_velx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          gf_vely(p.I) = 0.0;
          gf_velz(p.I) = 0.0;
          gf_press(p.I) = 1.0; // should add kinetic energy here
                               // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I) / (gf_rho(p.I) * (gamma - 1));
        });

  } else if (CCTK_EQUALS(initial_hydro, "shock tube")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            gf_rho(p.I) = 2.0;
            gf_velx(p.I) = 0.0;
            gf_vely(p.I) = 0.0;
            gf_velz(p.I) = 0.0;
            gf_press(p.I) = 2.0;
          } else {
            gf_rho(p.I) = 1.0;
            gf_velx(p.I) = 0.0;
            gf_vely(p.I) = 0.0;
            gf_velz(p.I) = 0.0;
            gf_press(p.I) = 1.0;
          }

          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I) / (gf_rho(p.I) * (gamma - 1));
        });

  } else if (CCTK_EQUALS(initial_hydro, "balsara1")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            gf_rho(p.I) = 1.0;
            gf_velx(p.I) = 0.0;
            gf_vely(p.I) = 0.0;
            gf_velz(p.I) = 0.0;
            gf_press(p.I) = 1.0;
          } else {
            gf_rho(p.I) = 0.125;
            gf_velx(p.I) = 0.0;
            gf_vely(p.I) = 0.0;
            gf_velz(p.I) = 0.0;
            gf_press(p.I) = 0.1;
          }

          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I) / (gf_rho(p.I) * (gamma - 1));
        });

  } else if (CCTK_EQUALS(initial_hydro, "spherical shock")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL r2 = pow2(p.x) + pow2(p.y) + pow2(p.z);
          if (r2 <= pow2(shock_radius)) {
            gf_rho(p.I) = 2.0;
            gf_velx(p.I) = 0.0;
            gf_vely(p.I) = 0.0;
            gf_velz(p.I) = 0.0;
            gf_press(p.I) = 2.0;
          } else {
            gf_rho(p.I) = 1.0;
            gf_velx(p.I) = 0.0;
            gf_vely(p.I) = 0.0;
            gf_velz(p.I) = 0.0;
            gf_press(p.I) = 1.0;
          }

          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I) / (gf_rho(p.I) * (gamma - 1));
        });

  } else {
    CCTK_ERROR("Internal error");
  }
}

} // namespace HydroInitial
