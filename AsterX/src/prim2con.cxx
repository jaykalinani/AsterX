#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "prim2con.hxx"

#include "eos.hxx"
#include "eos_idealgas.hxx"

namespace AsterX {
using namespace Loop;
using namespace EOSX;
using namespace std;

extern "C" void AsterX_Prim2Con_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Prim2Con_Initial;
  DECLARE_CCTK_PARAMETERS;

  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
       rgye(ye_min, ye_max);
  
  const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                   calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                   calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};

        if (entropy(p.I) != entropy(p.I) ||
            eps(p.I) != eps(p.I) ) {
          printf("Initial data is craaaaap");
          printf("rho = %16.8e \n",rho(p.I));
          printf("eps = %16.8e \n",eps(p.I));
          printf("press = %16.8e \n",press(p.I));
          printf("entropy = %16.8e \n",entropy(p.I));
          printf("velx = %16.8e \n",velx(p.I));
          printf("vely = %16.8e \n",vely(p.I));
          printf("velz = %16.8e \n",velz(p.I));
          printf("Bvecx = %16.8e \n",Bvecx(p.I));
          printf("Bvecy = %16.8e \n",Bvecy(p.I));
          printf("Bvecz = %16.8e \n",Bvecz(p.I));
          assert(0);
        }

        prim pv;
        pv.rho = rho(p.I);
        pv.vel(0) = velx(p.I);
        pv.vel(1) = vely(p.I);
        pv.vel(2) = velz(p.I);
        pv.eps = eps(p.I);
        pv.press = press(p.I);
        pv.Bvec(0) = Bvecx(p.I);
        pv.Bvec(1) = Bvecy(p.I);
        pv.Bvec(2) = Bvecz(p.I);

        // Lower velocity
        vec<CCTK_REAL, 3> v_low = calc_contraction(g, pv.vel);

        // Check validity of primitives ----------
        // Copy of code from c2p::prims_floors_and_ceilings
        
        if(true){

          const CCTK_REAL w_lim = sqrt(1.0 + vw_lim * vw_lim);
          const CCTK_REAL v_lim = vw_lim / w_lim;
      
          // ----------
          // Floor and ceiling for rho and velocity
          // Keeps pressure the same and changes eps
          // ----------
      
          // check if computed velocities are within the specified limit
          CCTK_REAL vsq_Sol = calc_contraction(v_low, pv.vel);
          CCTK_REAL sol_v = sqrt(vsq_Sol);
          if (sol_v > v_lim) {

            pv.vel *= v_lim / sol_v;
            v_low = v_low * v_lim / sol_v;
            velx(p.I) = pv.vel(0);
            vely(p.I) = pv.vel(1);
            velz(p.I) = pv.vel(2);
     
          }
      
          if (pv.rho > eos_th.rgrho.max) {
      
            // remove mass
            pv.rho = eos_th.rgrho.max;
            rho(p.I) = pv.rho;

            pv.eps = eos_th.eps_from_valid_rho_press_ye(eos_th.rgrho.max, pv.press, 0.5);
            eps(p.I) = pv.eps;

            entropy(p.I) =
                eos_th.kappa_from_valid_rho_eps_ye(eos_th.rgrho.max, pv.eps, 0.5);

          }

          if (pv.rho < rho_abs_min * (1 + atmo_tol)) {
      
            // add mass
            pv.rho = rho_abs_min;
            rho(p.I) = pv.rho;

            pv.eps = eos_th.eps_from_valid_rho_press_ye(rho_abs_min, pv.press, 0.5);
            eps(p.I) = pv.eps;

            entropy(p.I) =
                eos_th.kappa_from_valid_rho_eps_ye(rho_abs_min, pv.eps, 0.5);

          }
      
          // ----------
          // Floor and ceiling for eps
          // Keeps rho the same and changes press
          // ----------
      
          // check the validity of the computed eps
          auto rgeps = eos_th.range_eps_from_valid_rho_ye(pv.rho, 0.5);
          if (pv.eps > rgeps.max) {

            pv.eps = rgeps.max;
            eps(p.I) = pv.eps;

            pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, 0.5);
            press(p.I) = pv.press;

            entropy(p.I) = eos_th.kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, 0.5);  

          } else if (pv.eps < rgeps.min) {

            pv.eps = rgeps.min;
            eps(p.I) = pv.eps;

            pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, 0.5);
            press(p.I) = pv.press;

            entropy(p.I) = eos_th.kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, 0.5); 

          }

        }

        CCTK_REAL v2_new = calc_contraction(v_low,pv.vel); 
        if (pv.vel(0)>=1.0 || pv.vel(1)>=1.0 || pv.vel(2)>=1.0 || v2_new >= 1.0) {
          printf("Velocities are superluminal from ID generation. Abort.");
          assert(0);
        }
  
        // ---------- End of validity check

        cons cv;
        prim2con(g, pv, cv);

        dens(p.I) = cv.dens;
        momx(p.I) = cv.mom(0);
        momy(p.I) = cv.mom(1);
        momz(p.I) = cv.mom(2);
        tau(p.I) = cv.tau;
        dBx(p.I) = cv.dBvec(0);
        dBy(p.I) = cv.dBvec(1);
        dBz(p.I) = cv.dBvec(2);

        if (cv.tau!=cv.tau) {
          printf("NaN in cvtau - 1");
          assert(0);
        }

        if (max(max(abs(cv.dBvec(0)),abs(cv.dBvec(1))),abs(cv.dBvec(2))) > 10.0) {
          printf("B field too large.");
          printf("dBx = %16.8e \n",cv.dBvec(0));
          printf("dBy = %16.8e \n",cv.dBvec(1));
          printf("dBz = %16.8e \n",cv.dBvec(2));
          assert(0);
        }

       if (cv.dBvec(0)!=cv.dBvec(0)) {
          printf("NaN in bx");
          assert(0);
        }

        if (cv.dBvec(1)!=cv.dBvec(1)) {
          printf("NaN in by");
          assert(0);
        }

        if (cv.dBvec(2)!=cv.dBvec(2)) {
          printf("NaN in bz");
          assert(0);
        }
 

        DEnt(p.I) = entropy(p.I)*cv.dens;

        if (DEnt(p.I)!=DEnt(p.I)) {
          printf("NaN in DEnt in prim2con initial - chabamamba");
          assert(0);
        }

        if (tau(p.I)!=tau(p.I)) {
          printf("NaN in cvtau - 2");
          assert(0);
        }

        saved_rho(p.I) = pv.rho;
        saved_velx(p.I) = pv.vel(0);
        saved_vely(p.I) = pv.vel(1);
        saved_velz(p.I) = pv.vel(2);
        saved_eps(p.I) = pv.eps;

	const vec<CCTK_REAL, 3> v_up{pv.vel(0),pv.vel(1),pv.vel(2)};
        //const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
        CCTK_REAL wlor = calc_wlorentz(v_low, v_up);

	zvec_x(p.I) = wlor * pv.vel(0);
	zvec_y(p.I) = wlor * pv.vel(1);
	zvec_z(p.I) = wlor * pv.vel(2);

	svec_x(p.I) = (pv.rho+pv.rho*pv.eps+pv.press)*wlor*wlor*pv.vel(0);
	svec_y(p.I) = (pv.rho+pv.rho*pv.eps+pv.press)*wlor*wlor*pv.vel(1);
	svec_z(p.I) = (pv.rho+pv.rho*pv.eps+pv.press)*wlor*wlor*pv.vel(2);

      });

  if (zero_psi) {
    /* Initilaize Psi to 0.0 */
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                              CCTK_ATTRIBUTE_ALWAYS_INLINE { Psi(p.I) = 0.0; });
  }
}

} // namespace AsterX
