#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdio>
#include <cstdbool>
#include <cmath> 

#include "FM_disk_implementation.hxx"

extern "C" void FishboneMoncrief_ET_GRHD_initial(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTSX_FishboneMoncrief_ET_GRHD_initial;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Fishbone-Moncrief Disk Initial data");
  CCTK_VINFO("Using input parameters of\n a = %e,\n M = %e,\nr_in = %e,\nr_at_max_density = %e\nkappa = %e\ngamma = %e",a,M,r_in,r_at_max_density,kappa,gamma);

  const CCTK_REAL xcoord_init = r_at_max_density;
  const CCTK_REAL ycoord_init = 0.0;
  const CCTK_REAL zcoord_init = 0.0;
  const CCTK_REAL rr_init = sqrt(xcoord_init*xcoord_init+ycoord_init*ycoord_init+zcoord_init*zcoord_init);

  // First compute maximum pressure and density
  const CCTK_REAL hm1_init = FMdisk::GRHD_hm1(xcoord_init,ycoord_init,zcoord_init);
  const CCTK_REAL rho_max = pow( hm1_init * (gamma-1.0) / (kappa*gamma), 1.0/(gamma-1.0) );
  const CCTK_REAL P_max   = kappa * pow(rho_max, gamma);

  // We enforce units such that rho_max = 1.0; if these units are not obeyed, then
  //    we error out. If we did not error out, then the value of kappa used in all
  //    EOS routines would need to be changed, and generally these appear as
  //    read-only parameters.
  if(fabs(P_max/rho_max - kappa) > 1e-8) {
    printf("Error: To ensure that P = kappa*rho^Gamma, where rho_max = 1.0,\n");
    printf("       you must set (in your parfile) the polytropic constant kappa = P_max/rho_max = %.15e\n\n",P_max/rho_max);
    printf(" Needed values for kappa, for common values of Gamma:\n");
    printf(" For Gamma =4/3, use kappa=K_initial=K_poly = 4.249572342020724e-03 to ensure rho_max = 1.0\n");
    printf(" For Gamma =5/3, use kappa=K_initial=K_poly = 6.799315747233158e-03 to ensure rho_max = 1.0\n");
    printf(" For Gamma = 2,  use kappa=K_initial=K_poly = 8.499144684041449e-03 to ensure rho_max = 1.0\n");
    exit(1);
  }

  // Setup metric, loop over vertices
  grid.loop_int<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        CCTK_REAL xcoord = p.x;
        CCTK_REAL ycoord = p.y;
        CCTK_REAL zcoord = p.z;
        CCTK_REAL rr = sqrt(xcoord*xcoord+ycoord*ycoord+zcoord*zcoord);

        CCTK_REAL alp_L{0.};
        CCTK_REAL betaU0_L{0.};
        CCTK_REAL betaU1_L{0.};
        CCTK_REAL betaU2_L{0.};
        CCTK_REAL gDD00_L{0.};
        CCTK_REAL gDD01_L{0.};
        CCTK_REAL gDD02_L{0.};
        CCTK_REAL gDD11_L{0.};
        CCTK_REAL gDD12_L{0.};
        CCTK_REAL gDD22_L{0.};
        CCTK_REAL kDD00_L{0.};
        CCTK_REAL kDD01_L{0.};
        CCTK_REAL kDD02_L{0.};
        CCTK_REAL kDD11_L{0.};
        CCTK_REAL kDD12_L{0.};
        CCTK_REAL kDD22_L{0.};

        FMdisk::KerrSchild(xcoord,ycoord,zcoord,
                           alp_L,betaU0_L,betaU1_L,betaU2_L,
                           gDD00_L,gDD01_L,gDD02_L,gDD11_L,gDD12_L,gDD22_L,
                           kDD00_L,kDD01_L,kDD02_L,kDD11_L,kDD12_L,kDD22_L);

        alp(p.I) = alp_L;
        betax(p.I) = betaU0_L;
        betay(p.I) = betaU1_L;
        betaz(p.I) = betaU2_L;

        dtalp(p.I) = 0.0;
        dtbetax(p.I) = 0.0;
        dtbetay(p.I) = 0.0;
        dtbetaz(p.I) = 0.0;

        gxx(p.I) = gDD00_L;
        gyy(p.I) = gDD11_L;
        gzz(p.I) = gDD22_L;

        gxy(p.I) = gDD01_L;
        gxz(p.I) = gDD02_L;
        gyz(p.I) = gDD12_L;

        kxx(p.I) = kDD00_L;
        kyy(p.I) = kDD11_L;
        kzz(p.I) = kDD22_L;

        kxy(p.I) = kDD01_L;
        kxz(p.I) = kDD02_L;
        kyz(p.I) = kDD12_L;

  });

  // Setup hydro, loop over cell centers
  grid.loop_int<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        CCTK_REAL xcoord = p.x;
        CCTK_REAL ycoord = p.y;
        CCTK_REAL zcoord = p.z;
        CCTK_REAL rr = sqrt(xcoord*xcoord+ycoord*ycoord+zcoord*zcoord);

        bool set_to_atmosphere=false;

        if(rr > r_in) {

          CCTK_REAL hm1 = FMdisk::GRHD_hm1(xcoord,ycoord,zcoord);

          if(hm1 > 0) {

            rho(p.I) = pow( hm1 * (gamma-1.0) / (kappa*gamma), 1.0/(gamma-1.0) ) / rho_max;
            press(p.I) = kappa*pow(rho(p.I), gamma);
            // P = (\Gamma - 1) rho epsilon
            eps(p.I) = press(p.I) / (rho(p.I) * (gamma - 1.0));

            CCTK_REAL velU0_L{0.};
            CCTK_REAL velU1_L{0.};
            CCTK_REAL velU2_L{0.};

            FMdisk::GRHD_velocities(xcoord,ycoord,zcoord,
                                    velU0_L,velU1_L,velU2_L);

            velx(p.I) = velU0_L;
            vely(p.I) = velU1_L;
            velz(p.I) = velU2_L;

          } else {

            set_to_atmosphere=true;

          }

        } else {

          set_to_atmosphere=true;

        }
        // Outside the disk? Set to atmosphere all hydrodynamic variables!
        if(set_to_atmosphere) {
          // Choose an atmosphere such that
          //   rho =       1e-5 * r^(-3/2), and
          //   P   = k rho^gamma
          // Add 1e-100 or 1e-300 to rr or rho to avoid divisions by zero.
          rho(p.I) = 1e-5 * pow(rr + 1e-100,-3.0/2.0);
          press(p.I) = kappa*pow(rho(p.I), gamma);
          eps(p.I) = press(p.I) / ((rho(p.I) + 1e-300) * (gamma - 1.0));
          velx(p.I) = 0.0;
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
        }
      });

  /*
  CCTK_INT final_idx = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-1,cctk_lsh[1]-1,cctk_lsh[2]-1);
  CCTK_VINFO("=====   OUTPUTS   =====");
  CCTK_VINFO("betai: %e %e %e \ngij: %e %e %e %e %e %e \nKij: %e %e %e %e %e %e\nalp: %e\n",
              betax[final_idx],betay[final_idx],betaz[final_idx],
              gxx[final_idx],gxy[final_idx],gxz[final_idx],gyy[final_idx],gyz[final_idx],gzz[final_idx],
              kxx[final_idx],kxy[final_idx],kxz[final_idx],kyy[final_idx],kyz[final_idx],kzz[final_idx],
              alp[final_idx]);
  CCTK_VINFO("rho: %.15e\nPressure: %.15e\nvx: %.15e\nvy: %.15e\nvz: %.15e",
              rho[final_idx],press[final_idx],
              velx[final_idx],vely[final_idx],velz[final_idx]);
  */
}

extern "C" void FishboneMoncrief_ET_GRHD_initial__perturb_pressure(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTSX_FishboneMoncrief_ET_GRHD_initial__perturb_pressure;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

    CCTK_REAL xcoord = p.x;
    CCTK_REAL ycoord = p.y;
    CCTK_REAL zcoord = p.z;
    CCTK_REAL rr = sqrt(xcoord*xcoord+ycoord*ycoord+zcoord*zcoord);

    CCTK_REAL hm1 = 0.;

    if(rr > r_in) {

	hm1 = FMdisk::GRHD_hm1(xcoord,ycoord,zcoord);

        if(hm1 > 0) {

          CCTK_REAL press_L = press(p.I);
          CCTK_REAL eps_L = eps(p.I);
          CCTK_REAL rho_L = rho(p.I);

          FMdisk::GRHD_perturb_pressure(press_L, eps_L, rho_L);

          press(p.I) = press_L;
          eps(p.I) = eps_L;
          rho(p.I) = rho_L;

        }
    }
  });
}

extern "C" void FishboneMoncrief_Set_Spacetime(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTSX_FishboneMoncrief_Set_Spacetime;
  DECLARE_CCTK_PARAMETERS;

  // Setup metric, loop over vertices
  grid.loop_int<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

        CCTK_REAL xcoord = p.x;
        CCTK_REAL ycoord = p.y;
        CCTK_REAL zcoord = p.z;
        CCTK_REAL rr = sqrt(xcoord*xcoord+ycoord*ycoord+zcoord*zcoord);

        CCTK_REAL alp_L{0.};
        CCTK_REAL betaU0_L{0.};
        CCTK_REAL betaU1_L{0.};
        CCTK_REAL betaU2_L{0.};
        CCTK_REAL gDD00_L{0.};
        CCTK_REAL gDD01_L{0.};
        CCTK_REAL gDD02_L{0.};
        CCTK_REAL gDD11_L{0.};
        CCTK_REAL gDD12_L{0.};
        CCTK_REAL gDD22_L{0.};
        CCTK_REAL kDD00_L{0.};
        CCTK_REAL kDD01_L{0.};
        CCTK_REAL kDD02_L{0.};
        CCTK_REAL kDD11_L{0.};
        CCTK_REAL kDD12_L{0.};
        CCTK_REAL kDD22_L{0.};

        FMdisk::KerrSchild(xcoord,ycoord,zcoord,
                           alp_L,betaU0_L,betaU1_L,betaU2_L,
                           gDD00_L,gDD01_L,gDD02_L,gDD11_L,gDD12_L,gDD22_L,
                           kDD00_L,kDD01_L,kDD02_L,kDD11_L,kDD12_L,kDD22_L);

        alp(p.I) = alp_L;
        betax(p.I) = betaU0_L;
        betay(p.I) = betaU1_L;
        betaz(p.I) = betaU2_L;

        dtalp(p.I) = 0.0;
        dtbetax(p.I) = 0.0;
        dtbetay(p.I) = 0.0;
        dtbetaz(p.I) = 0.0;

        gxx(p.I) = gDD00_L;
        gyy(p.I) = gDD11_L;
        gzz(p.I) = gDD22_L;

        gxy(p.I) = gDD01_L;
        gxz(p.I) = gDD02_L;
        gyz(p.I) = gDD12_L;

        kxx(p.I) = kDD00_L;
        kyy(p.I) = kDD11_L;
        kzz(p.I) = kDD22_L;

        kxy(p.I) = kDD01_L;
        kxz(p.I) = kDD02_L;
        kyz(p.I) = kDD12_L;

  });
}