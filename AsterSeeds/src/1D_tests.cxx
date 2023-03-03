#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <seeds_utils.hxx>

#include <eos.hxx>
#include <eos_idealgas.hxx>

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace EOSX;

extern "C" void Tests1D_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Tests1D_Initialize;
  DECLARE_CCTK_PARAMETERS;


  // For all the tests, the initial data EOS is ideal gas
  // Constructing the IG EOS object
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);
  const CCTK_REAL dummy_ye = 0.5;

  if (CCTK_EQUALS(test_case, "equilibrium")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0;
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I), dummy_ye);
          // Bvecx(p.I) = 0.0;
          // Bvecy(p.I) = 0.0;
          // Bvecz(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(test_case, "sound wave")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0; // should add kinetic energy here
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I), dummy_ye);
          // Bvecx(p.I) = 0.0;
          // Bvecy(p.I) = 0.0;
          // Bvecz(p.I) = 0.0;
        });

  } else if (CCTK_EQUALS(test_case, "shock tube")) {
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
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I), dummy_ye);
          // Bvecx(p.I) = 0.0;
          // Bvecy(p.I) = 0.0;
          // Bvecz(p.I) = 0.0;
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

  } else if (CCTK_EQUALS(test_case, "Balsara1")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 1.0;
            // Bvecx(p.I) = 0.5;
            // Bvecy(p.I) = 1.0;
            // Bvecz(p.I) = 0.0;

          } else {
            rho(p.I) = 0.125;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 0.1;
            // Bvecx(p.I) = 0.5;
            // Bvecy(p.I) = -1.0;
            // Bvecz(p.I) = 0.0;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I), dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(grid.nghostzones,
                                  [=] CCTK_DEVICE(const PointDesc &p)
                                      CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                        if (p.x <= 0.0) {
                                          Avec_x(p.I) = 1.0 * (p.z);
                                        } else {
                                          Avec_x(p.I) = -1.0 * (p.z);
                                        }
                                      });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.5 * (p.y); });

  } else if (CCTK_EQUALS(test_case, "Balsara2")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 30.0;
            // Bvecx(p.I) = 5.0;
            // Bvecy(p.I) = 6.0;
            // Bvecz(p.I) = 6.0;
          } else {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 1.0;
            // Bvecx(p.I) = 5.0;
            // Bvecy(p.I) = 0.7;
            // Bvecz(p.I) = 0.7;
          }
	  eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I), dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            Avec_x(p.I) = 6.0 * (p.z) - 6.0 * (p.y);
          } else {
            Avec_x(p.I) = 0.7 * (p.z) - 0.7 * (p.y);
          }
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 5.0 * (p.y); });

  } else if (CCTK_EQUALS(test_case, "Balsara3")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 1000.0;
            // Bvecx(p.I) = 10.0;
            // Bvecy(p.I) = 7.0;
            // Bvecz(p.I) = 7.0;

          } else {
            rho(p.I) = 1.0;
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 0.1;
            // Bvecx(p.I) = 10.0;
            // Bvecy(p.I) = 0.7;
            // Bvecz(p.I) = 0.7;
          }
	  eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I), dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            Avec_x(p.I) = 7.0 * (p.z) - 7.0 * (p.y);
          } else {
            Avec_x(p.I) = 0.7 * (p.z) - 0.7 * (p.y);
          }
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 10.0 * (p.y); });

  } else if (CCTK_EQUALS(test_case, "Balsara4")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 1.0;
            velx(p.I) = 0.999;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 0.1;
            // Bvecx(p.I) = 10.0;
            // Bvecy(p.I) = 7.0;
            // Bvecz(p.I) = 7.0;

          } else {
            rho(p.I) = 1.0;
            velx(p.I) = -0.999;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            press(p.I) = 0.1;
            // Bvecx(p.I) = 10.0;
            // Bvecy(p.I) = -7.0;
            // Bvecz(p.I) = -7.0;
          }
	  eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I), dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            Avec_x(p.I) = 7.0 * (p.z) - 7.0 * (p.y);
          } else {
            Avec_x(p.I) = -7.0 * (p.z) + 7.0 * (p.y);
          }
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 10.0 * (p.y); });

  } else if (CCTK_EQUALS(test_case, "Balsara5")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = 1.08;
            velx(p.I) = 0.4;
            vely(p.I) = 0.3;
            velz(p.I) = 0.2;
            press(p.I) = 0.95;
            // Bvecx(p.I) = 2.0;
            // Bvecy(p.I) = 0.3;
            // Bvecz(p.I) = 0.3;

          } else {
            rho(p.I) = 1.0;
            velx(p.I) = -0.45;
            vely(p.I) = -0.2;
            velz(p.I) = 0.2;
            press(p.I) = 1.0;
            // Bvecx(p.I) = 2.0;
            // Bvecy(p.I) = -0.7;
            // Bvecz(p.I) = 0.5;
          }
	  eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I), dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            Avec_x(p.I) = 0.3 * (p.z) - 0.3 * (p.y);
          } else {
            Avec_x(p.I) = -0.7 * (p.z) - 0.5 * (p.y);
          }
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 2.0 * (p.y); });

  } else {
    CCTK_ERROR("Test case not defined");
  }
}

} // namespace AsterSeeds
