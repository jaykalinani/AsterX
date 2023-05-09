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
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
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

  } else if (CCTK_EQUALS(test_case, "sound wave")) {
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0; // should add kinetic energy here
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
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
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
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
   
    CCTK_REAL rhol = 1.0;
    CCTK_REAL vxl = 0.0;
    CCTK_REAL vyl = 0.0;
    CCTK_REAL vzl = 0.0;
    CCTK_REAL pressl = 1.0;
    CCTK_REAL Bxl = 0.5;
    CCTK_REAL Byl = 1.0;
    CCTK_REAL Bzl = 0.0;

    CCTK_REAL rhor = 0.125;
    CCTK_REAL vxr = 0.0;
    CCTK_REAL vyr = 0.0;
    CCTK_REAL vzr = 0.0;
    CCTK_REAL pressr = 0.1;
    CCTK_REAL Bxr = 0.5;
    CCTK_REAL Byr = -1.0;
    CCTK_REAL Bzr = 0.0;    
	 
    if (CCTK_EQUALS(shock_dir, "x")) {
   
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.x <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = vxl;
            vely(p.I) = vyl;
            velz(p.I) = vzl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = vxr;
            vely(p.I) = vyr;
            velz(p.I) = vzr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.x <= 0.0) { 
	       Avec_x(p.I) = Byl * (p.z) - Bzl * (p.y);  
	    } else { 
	       Avec_x(p.I) = Byr * (p.z) - Bzr * (p.y); } });


    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = 0.0; });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.x <= 0.0) {
               Avec_z(p.I) = Bxl * (p.y);
            } else {
               Avec_z(p.I) = Bxr * (p.y); } });
    
    } else if (CCTK_EQUALS(shock_dir, "y")) {

    //x-->z, y-->x, z-->y
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.y <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = vzl;
            vely(p.I) = vxl;
            velz(p.I) = vyl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = vzr;
            vely(p.I) = vxr;
            velz(p.I) = vyr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.y <= 0.0) {
               Avec_z(p.I) = Byl * (p.z) - Bzl * (p.y);
            } else {
               Avec_z(p.I) = Byr * (p.z) - Bzr * (p.y); } });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_x(p.I) = 0.0; });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.y <= 0.0) {
               Avec_y(p.I) = Bxl * (p.y);
            } else {
               Avec_y(p.I) = Bxr * (p.y); } });

    } else if (CCTK_EQUALS(shock_dir, "z")) {
    
    //x-->y, y-->z, z-->x
    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          if (p.z <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = vyl;
            vely(p.I) = vzl;
            velz(p.I) = vxl;
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = vyr;
            vely(p.I) = vzr;
            velz(p.I) = vxr;
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all_device<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.z <= 0.0) {
               Avec_y(p.I) = Byl * (p.z) - Bzl * (p.y);
            } else {
               Avec_y(p.I) = Byr * (p.z) - Bzr * (p.y); } });

    grid.loop_all_device<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = 0.0; });

    grid.loop_all_device<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
            if (p.z <= 0.0) {
               Avec_x(p.I) = Bxl * (p.y);
            } else {
               Avec_x(p.I) = Bxr * (p.y); } });

    } else { CCTK_ERROR("Shock direction case not defined"); }
    

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
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
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
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
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
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
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
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
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
