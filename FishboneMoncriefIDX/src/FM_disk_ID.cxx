#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdio>
#include <cstdbool>
#include <cmath> 

#include "FM_disk_implementation.hxx"

// Alias for "vel" vector gridfunction:
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])

extern "C" void FishboneMoncrief_ET_GRHD_initial(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTSX_FishboneMoncrief_ET_GRHD_initial;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Fishbone-Moncrief Disk Initial data");
  CCTK_VINFO("Using input parameters of\n a = %e,\n M = %e,\nr_in = %e,\nr_at_max_density = %e\nkappa = %e\ngamma = %e",a,M,r_in,r_at_max_density,kappa,gamma);

  // First compute maximum pressure and density
  CCTK_REAL hm1 = FMdisk::GRHD_hm1();
  CCTK_REAL xcoord = r_at_max_density;
  CCTK_REAL ycoord = 0.0;
  CCTK_REAL zcoord = 0.0;
  CCTK_REAL rho_max = pow( hm1 * (gamma-1.0) / (kappa*gamma), 1.0/(gamma-1.0) );
  CCTK_REAL P_max   = kappa * pow(rho_max, gamma);

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

  grid.loop_int<1, 1, 1>(
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

        alp_cell(p.I) = alp_L;
        betax_cell(p.I) = betaU0_L;
        betay_cell(p.I) = betaU1_L;
        betaz_cell(p.I) = betaU2_L;

        dtalp_cell(p.I) = 0.0;
        dtbetax_cell(p.I) = 0.0;
        dtbetay_cell(p.I) = 0.0;
        dtbetaz_cell(p.I) = 0.0;

        gxx_cell(p.I) = gDD00_L;
        gyy_cell(p.I) = gDD11_L;
        gzz_cell(p.I) = gDD22_L;

        gxy_cell(p.I) = gDD01_L;
        gxz_cell(p.I) = gDD02_L;
        gyz_cell(p.I) = gDD12_L;

        kxx_cell(p.I) = kDD00_L;
        kyy_cell(p.I) = kDD11_L;
        kzz_cell(p.I) = kDD22_L;

        kxy_cell(p.I) = kDD01_L;
        kxz_cell(p.I) = kDD02_L;
        kyz_cell(p.I) = kDD12_L;

// HERE

        bool set_to_atmosphere=false;
        if(rr > r_in) {
          if(hm1 > 0) {

            rho[idx] = pow( hm1 * (gamma-1.0) / (kappa*gamma), 1.0/(gamma-1.0) ) / rho_max;
            press[idx] = kappa*pow(rho[idx], gamma);
            // P = (\Gamma - 1) rho epsilon
            eps[idx] = press[idx] / (rho[idx] * (gamma - 1.0));
            FishboneMoncrief_FMdisk_GRHD_velocities(cctkGH,cctk_lsh,
                                                    i,j,k,
                                                    x,y,z,
                                                    velx,vely,velz);
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
          rho[idx] = 1e-5 * pow(rr + 1e-100,-3.0/2.0);
          press[idx] = kappa*pow(rho[idx], gamma);
          eps[idx] = press[idx] / ((rho[idx] + 1e-300) * (gamma - 1.0));
          w_lorentz[idx] = 1.0;
          velx[idx] = 0.0;
          vely[idx] = 0.0;
          velz[idx] = 0.0;
        }
      }

  CCTK_INT final_idx = CCTK_GFINDEX3D(cctkGH,cctk_lsh[0]-1,cctk_lsh[1]-1,cctk_lsh[2]-1);
  CCTK_VINFO("=====   OUTPUTS   =====");
  CCTK_VINFO("betai: %e %e %e \ngij: %e %e %e %e %e %e \nKij: %e %e %e %e %e %e\nalp: %e\n",betax[final_idx],betay[final_idx],betaz[final_idx],gxx[final_idx],gxy[final_idx],gxz[final_idx],gyy[final_idx],gyz[final_idx],gzz[final_idx],kxx[final_idx],kxy[final_idx],kxz[final_idx],kyy[final_idx],kyz[final_idx],kzz[final_idx],alp[final_idx]);
  CCTK_VINFO("rho: %.15e\nPressure: %.15e\nvx: %.15e\nvy: %.15e\nvz: %.15e",rho[final_idx],press[final_idx],velx[final_idx],vely[final_idx],velz[final_idx]);
}
