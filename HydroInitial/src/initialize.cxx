#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace HydroInitial {
using namespace std;
using namespace Loop;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
}

extern "C" void HydroInitial_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_HydroInitial_Initialize;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(initial_hydro, "equilibrium")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0;
          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          eps(p.I) = press(p.I) / (rho(p.I) * (gamma - 1));
          //Bvecx(p.I) = 0.0;
          //Bvecy(p.I) = 0.0;
          //Bvecz(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(initial_hydro, "sound wave")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0; // should add kinetic energy here
                               // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          eps(p.I) = press(p.I) / (rho(p.I) * (gamma - 1));
          //Bvecx(p.I) = 0.0;
          //Bvecy(p.I) = 0.0;
          //Bvecz(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(initial_hydro, "shock tube")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 2.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 2.0;
          } else {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 1.0;
          }

          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          eps(p.I) = press(p.I) / (rho(p.I) * (gamma - 1));
          //Bvecx(p.I) = 0.0;
          //Bvecy(p.I) = 0.0;
          //Bvecz(p.I) = 0.0;
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
	      Avec_x(p.I) = 0.0;
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          Avec_y(p.I) = 0.0;
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          Avec_z(p.I) = 0.0;
        });


  } else if (CCTK_EQUALS(initial_hydro, "Balsara1")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 1.0;
            //Bvecx(p.I) = 0.5;
            //Bvecy(p.I) = 1.0;
            //Bvecz(p.I) = 0.0;

          } else {
            rho(p.I) = 0.125;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 0.1;
            //Bvecx(p.I) = 0.5;
            //Bvecy(p.I) = -1.0;
            //Bvecz(p.I) = 0.0;
          }
          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          eps(p.I) = press(p.I) / (rho(p.I) * (gamma - 1));
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if(p.x <= 0.0) {
            Avec_x(p.I) = 1.0 * (p.z);
          } else {
            Avec_x(p.I) = -1.0 * (p.z);
          }
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          Avec_y(p.I) = 0.0;
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          Avec_z(p.I) = 0.5 * (p.y);
        });

  } else if (CCTK_EQUALS(initial_hydro, "spherical shock")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const CCTK_REAL r2 = pow2(p.x) + pow2(p.y) + pow2(p.z);
          if (r2 <= pow2(shock_radius)) {
            rho(p.I) = 2.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 2.0;
          } else {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 1.0;
          }

          // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          eps(p.I) = press(p.I) / (rho(p.I) * (gamma - 1));
          //Bvecx(p.I) = 0.0;
          //Bvecy(p.I) = 0.0;
          //Bvecz(p.I) = 0.0;
        });

  } else {
    CCTK_ERROR("Internal error");
  }
}

} // namespace HydroInitial
