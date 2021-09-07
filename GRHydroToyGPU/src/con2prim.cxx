#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "reprimand/eos_thermal.h"  // The EOS framework
#include "reprimand/eos_idealgas.h"
#include "reprimand/con2prim_imhd.h" // The con2prim framework
using namespace EOS_Toolkit;

namespace GRHydroToyGPU {
    using namespace std;
    using namespace Loop;

// TODO: remove this if not needed
/*    namespace {
        template <typename T>
        inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
            return x*x;
        }
    }
*/

    /***************************************************************************
        What follows has been inspired by the content of /example/minimal.cc,
        which can be found inside the RePrimAnd library:
        https://github.com/wokast/RePrimAnd
    ****************************************************************************/

    // RePrimAnd C2P
    extern "C" void GRHydroToyGPU_Con2Prim_RePrimAnd(CCTK_ARGUMENTS) {
        DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Con2Prim_RePrimAnd;
        DECLARE_CCTK_PARAMETERS;

        // CCTK_INT di, dj, dk;

        const GridDescBaseDevice grid(cctkGH);
        constexpr array<int, dim> cell_centred   = {1, 1, 1};
        constexpr array<int, dim> vertex_centred = {0, 0, 0};
        const GF3D2layout gf_layout_cell(cctkGH,   cell_centred);
        const GF3D2layout gf_layout_vertex(cctkGH, vertex_centred);

        const GF3D2<const CCTK_REAL> gf_gxx(gf_layout_vertex, gxx);
        const GF3D2<const CCTK_REAL> gf_gxy(gf_layout_vertex, gxy);
        const GF3D2<const CCTK_REAL> gf_gxz(gf_layout_vertex, gxz);
        const GF3D2<const CCTK_REAL> gf_gyy(gf_layout_vertex, gyy);
        const GF3D2<const CCTK_REAL> gf_gyz(gf_layout_vertex, gyz);
        const GF3D2<const CCTK_REAL> gf_gzz(gf_layout_vertex, gzz);

        // CCTK_REAL gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg;

        GF3D2<CCTK_REAL> gf_dens(gf_layout_cell,  dens);
        GF3D2<CCTK_REAL> gf_momx(gf_layout_cell,  momx);
        GF3D2<CCTK_REAL> gf_momy(gf_layout_cell,  momy);
        GF3D2<CCTK_REAL> gf_momz(gf_layout_cell,  momz);
        GF3D2<CCTK_REAL> gf_tau(gf_layout_cell,   tau);

        // FIXME: import primitives from HydroBase? And HydroBase::vel is a 3-vector!
        GF3D2<CCTK_REAL> gf_rho(gf_layout_cell,   rho);
        GF3D2<CCTK_REAL> gf_velx(gf_layout_cell,  velx);
        GF3D2<CCTK_REAL> gf_vely(gf_layout_cell,  vely);
        GF3D2<CCTK_REAL> gf_velz(gf_layout_cell,  velz);
        GF3D2<CCTK_REAL> gf_press(gf_layout_cell, press);
        GF3D2<CCTK_REAL> gf_eps(gf_layout_cell,   eps);

        // TODO: add magnetic fields and tabulated EOS


        // Setting up the EOS
        /*const CCTK_REAL max_eps   = 11.;  // TODO *maybe*: real_t in place of CCTK_REAL
        const CCTK_REAL max_rho   = 1.e6;
        const CCTK_REAL adiab_ind = 1.;*/
        //const auto eos = make_eos_idealgas(adiab_ind, max_eps, max_rho);
        const auto eos = make_eos_idealgas(gl_gamma, gl_max_eps, gl_max_rho);

        //const eos_thermal& eos = global_eos_thermal::get_eos();
        //const eos_barotr& eos_id = global_eos_cold::get_eos();

        // Setting up artificial atmosphere
        /*const CCTK_REAL atmo_rho = 1.e-20;
        const CCTK_REAL atmo_eps = 0.1;
        const CCTK_REAL atmo_ye  = 0.5;
        const CCTK_REAL atmo_cut = atmo_rho*1.01;*/
        //const CCTK_REAL atmo_p   = eos.at_rho_eps_ye(atmo_rho, atmo_eps, atmo_ye).press();
        const CCTK_REAL atmo_p = eos.at_rho_eps_ye(atmo_rho, atmo_eps, atmo_ye).press();

        atmosphere atmo{atmo_rho, atmo_eps, atmo_ye, atmo_p, atmo_cut};

        /*
        CCTK_REAL rho_atmo_cut = rho_abs_min*(1.0 + Spritz_atmo_tolerance);
        assert(eos.is_rho_valid(rho_abs_min));
        assert(eos.is_ye_valid(Spritz_hot_atmo_Y_e));
        //CCTK_REAL eps_atmo    = eos_id.at_rho_ye(rho_abs_min,Spritz_hot_atmo_Y_e).eps();
        CCTK_REAL eps_atmo    = eos_id.at_rho(rho_abs_min).eps();
        eps_atmo    = eos.range_eps(rho_abs_min,Spritz_hot_atmo_Y_e).limit_to(eps_atmo);
        CCTK_REAL p_atmo  = eos.at_rho_eps_ye(rho_abs_min, eps_atmo, Spritz_hot_atmo_Y_e).press();
        const atmosphere atmo(rho_abs_min, eps_atmo, Spritz_hot_atmo_Y_e, p_atmo, rho_atmo_cut); 
        */


        // Primitive recovery parameters
        /*
        const CCTK_REAL    rho_strict = 1.e-11;
        const CCTK_BOOLEAN ye_lenient = false;
        const CCTK_INT     max_iter   = 30;
        const CCTK_REAL    c2p_acc    = 1.e-8;
        const CCTK_REAL    max_b      = 10.;
        const CCTK_REAL    max_z      = 1.e3;
        */
  
        // Get a recovery function
        con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, max_z,
                           max_b, atmo, c2p_acc, max_iter);

        /*           
        // Initializing con2prim object for inside-BH-horizon evolution 
        const con2prim_mhd cv2pv_BH (eos, rho_strict_BH, y_e_linient_BH,
                              maximum_z_BH, maximum_b_BH, atmo, Spritz_mhd_tolf, Spritz_countmax);

        //to validate storage when *GetTemp_FromThermo!=0
        if (*GetTemp_FromThermo) 
        {
         assert(temperature);
         assert(entropy);
        }
        */



        // Loop over the entire grid (0 to n-1 points in each direction)
        grid.loop_all_device<1, 1, 1>(
            grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
                CCTK_ATTRIBUTE_ALWAYS_INLINE {

            // Build the objects needed (primitive variables and error report)
            prim_vars_mhd pv;
            con2prim_mhd::report rep;

	    //dummy variables
            CCTK_REAL w_lorentz = 1.0;
            CCTK_REAL Y_e = 0.5;
            CCTK_REAL Y_e_con = 0.5*gf_dens(p.I);
	    CCTK_REAL Bx = 0.0;
	    CCTK_REAL By = 0.0;
	    CCTK_REAL Bz = 0.0;
	    CCTK_REAL dBx = 0.0;
	    CCTK_REAL dBy = 0.0;
	    CCTK_REAL dBz = 0.0;
	    CCTK_REAL Ex = 0.0;
	    CCTK_REAL Ey = 0.0;
	    CCTK_REAL Ez = 0.0;
	    
            // Interpolate metric terms from vertices to center
            CCTK_REAL gxx_avg = 0.;
            CCTK_REAL gxy_avg = 0.;
            CCTK_REAL gxz_avg = 0.;
            CCTK_REAL gyy_avg = 0.;
            CCTK_REAL gyz_avg = 0.;
            CCTK_REAL gzz_avg = 0.;

            for(int dk = 0 ; dk < 2 ; ++dk)
                for(int dj = 0 ; dj < 2 ; ++dj)
                    for(int di = 0 ; di < 2 ; ++di) {
                        gxx_avg += gf_gxx(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gxy_avg += gf_gxy(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gxz_avg += gf_gxz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gyy_avg += gf_gyy(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gyz_avg += gf_gyz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                        gzz_avg += gf_gzz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                    }

            gxx_avg *= 0.125;
            gxy_avg *= 0.125;
            gxz_avg *= 0.125;
            gyy_avg *= 0.125;
            gyz_avg *= 0.125;
            gzz_avg *= 0.125;

            sm_metric3 g(sm_symt3l(gxx_avg, gxy_avg, gyy_avg, gxz_avg, gyz_avg, gzz_avg));


            /* Build the object containing the conservative variables and
               perform the con2prim                                             */
            cons_vars_mhd cv{gf_dens(p.I), gf_tau(p.I), Y_e_con, 
                            {gf_momx(p.I), gf_momy(p.I), gf_momz(p.I)},
                            {dBx, dBy, dBz}};
            cv2pv(pv, cv, g, rep);


            // Handle incorrectable errors
            if (rep.failed()) {
                CCTK_WARN(1, rep.debug_message().c_str());

/*		
                // ***** DEBUG LINES ***** 
                CCTK_WARN(1,  "rep.failed() = True, printing variables BEFORE setting to atmosphere");
                CCTK_VWARN(1, "dx = %f", CCTK_DELTA_SPACE(0));
                CCTK_VWARN(1, "cctk_levfac = %f", double(cctk_levfac[0]));
                CCTK_VWARN(1, "cctk_iteration = %i", cctk_iteration);
                CCTK_VWARN(1, "cctk_time = %f", cctk_time);
                CCTK_VWARN(1, "p.I = %26.16e", p.I);
                //CCTK_VWARN(1, "i, j, k (in C) = %d, %d, %d\n", i, j, k);
                CCTK_VWARN(1, "x[p.I], y[p.I], z[p.I] = %26.16e, %26.16e, %26.16e\n",
                               x[p.I], y[p.I], z[p.I]);
                CCTK_VWARN(1, "g.vol_elem = %26.16e", g.vol_elem);
                CCTK_VWARN(1, "dens[p.I] = %26.16e", gf_dens[p.I]);
                CCTK_VWARN(1, "tau[p.I] = %26.16e", gf_tau[p.I]);
                CCTK_VWARN(1, "momx[p.I] = %26.16e", gf_momx[p.I]);
                CCTK_VWARN(1, "momy[p.I] = %26.16e", gf_momy[p.I]);
                CCTK_VWARN(1, "momz[p.I] = %26.16e", gf_momz[p.I]);
                
                CCTK_VWARN(1, "dBx[p.I] = %26.16e\n",     dBx[ijk]);
                CCTK_VWARN(1, "dBy[p.I] = %26.16e\n",     dBy[ijk]);
                CCTK_VWARN(1, "dBz[p.I] = %26.16e\n",     dBz[ijk]);
                CCTK_VWARN(1, "Y_e_con[p.I] = %26.16e\n", Y_e_con[ijk]);
                

                
                cv.bcons(0) = dBx[ijk];
                cv.bcons(1) = dBy[ijk];
                cv.bcons(2) = dBz[ijk];
                pv.B = cv.bcons / g.vol_elem;
                
*/
		atmo.set(pv, cv, g);
            }

            // Write back primitive variables
            pv.scatter(gf_rho(p.I), gf_eps(p.I), Y_e, gf_press(p.I),
                       gf_velx(p.I), gf_vely(p.I), gf_velz(p.I), w_lorentz,
                       Ex, Ey, Ez,
                       Bx, By, Bz);

            // Computing temperature and entropy
            /*     if(*GetTemp_FromThermo!=0)
            {  
                 auto state = eos.at_rho_eps_ye(rho[ijk], eps[ijk], Y_e[ijk]);
                 assert(state);
                 temperature[ijk] = state.temp();
                 entropy[ijk]     = state.sentr();
            } */
             
            /* Write back conserved variables in case they have been adjusted
               or set to NaN                                                    */
            if (rep.adjust_cons)
                cv.scatter(gf_dens(p.I), gf_tau(p.I), Y_e_con, 
                           gf_momx(p.I), gf_momy(p.I), gf_momy(p.I),
                           dBx, dBy, dBz); 
        });
    }
} // namespace GRHydroToyGPU



