#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "reprimand/eos_thermal.h" // The EOS framework
#include "reprimand/eos_idealgas.h"
#include "reprimand/con2prim_imhd.h" // The con2prim framework
using namespace EOS_Toolkit;

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

/***************************************************************************
    What follows has been inspired by the content of /example/minimal.cc,
    which can be found inside the RePrimAnd library:
    https://github.com/wokast/RePrimAnd
****************************************************************************/

// RePrimAnd C2P
extern "C" void GRHydroToyGPU_Con2Prim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  constexpr auto DI = PointDesc::DI;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  constexpr array<int, dim> vertex_centred = {0, 0, 0};
  const GF3D2layout gf_layout_cell(cctkGH, cell_centred);
  const GF3D2layout gf_layout_vertex(cctkGH, vertex_centred);

  const GF3D2<const CCTK_REAL> gf_gxx(gf_layout_vertex, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(gf_layout_vertex, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(gf_layout_vertex, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(gf_layout_vertex, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(gf_layout_vertex, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(gf_layout_vertex, gzz);

  GF3D2<CCTK_REAL> gf_dens(gf_layout_cell, dens);
  GF3D2<CCTK_REAL> gf_momx(gf_layout_cell, momx);
  GF3D2<CCTK_REAL> gf_momy(gf_layout_cell, momy);
  GF3D2<CCTK_REAL> gf_momz(gf_layout_cell, momz);
  GF3D2<CCTK_REAL> gf_tau(gf_layout_cell, tau);

  GF3D2<CCTK_REAL> gf_rho(gf_layout_cell, rho);
  GF3D2<CCTK_REAL> gf_velx(gf_layout_cell, velx);
  GF3D2<CCTK_REAL> gf_vely(gf_layout_cell, vely);
  GF3D2<CCTK_REAL> gf_velz(gf_layout_cell, velz);
  GF3D2<CCTK_REAL> gf_press(gf_layout_cell, press);
  GF3D2<CCTK_REAL> gf_eps(gf_layout_cell, eps);

  // TODO: add magnetic fields and tabulated EOS
  const CCTK_REAL adiab_ind = 1.0 / (gamma - 1);
  const auto eos = make_eos_idealgas(adiab_ind, max_eps, max_rho);
  const CCTK_REAL atmo_p =
      eos.at_rho_eps_ye(atmo_rho, atmo_eps, atmo_ye).press();
  atmosphere atmo{atmo_rho, atmo_eps, atmo_ye, atmo_p, atmo_cut};

  // Get a recovery function
  con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, max_z, max_b, atmo, c2p_acc,
                     max_iter);

  // Loop over the entire grid (0 to n-1 points in each direction)
  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Build the objects needed (primitive variables and error report)
        prim_vars_mhd pv;
        con2prim_mhd::report rep;

        // dummy variables
        CCTK_REAL w_lorentz = 1.0;
        CCTK_REAL Y_e = 0.5;
        CCTK_REAL Y_e_con = 0.5 * gf_dens(p.I);
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

        for (int dk = 0; dk < 2; ++dk)
          for (int dj = 0; dj < 2; ++dj)
            for (int di = 0; di < 2; ++di) {
              gxx_avg += gf_gxx(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              gxy_avg += gf_gxy(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              gxz_avg += gf_gxz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              gyy_avg += gf_gyy(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              gyz_avg += gf_gyz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              gzz_avg += gf_gzz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
            }

        gxx_avg *= 0.125;
        gxy_avg *= 0.125;
        gxz_avg *= 0.125;
        gyy_avg *= 0.125;
        gyz_avg *= 0.125;
        gzz_avg *= 0.125;

        sm_metric3 g(
            sm_symt3l(gxx_avg, gxy_avg, gyy_avg, gxz_avg, gyz_avg, gzz_avg));
        /* Build the object containing the conservative variables and
           perform the con2prim                                             */
        cons_vars_mhd cv{gf_dens(p.I),
                         gf_tau(p.I),
                         Y_e_con,
                         {gf_momx(p.I), gf_momy(p.I), gf_momz(p.I)},
                         {dBx, dBy, dBz}};
        cv2pv(pv, cv, g, rep);

        // Handle incorrectable errors
        if (rep.failed()) {
          CCTK_WARN(1, rep.debug_message().c_str());
          atmo.set(pv, cv, g);
        }

        // Write back primitive variables
        pv.scatter(gf_rho(p.I), gf_eps(p.I), Y_e, gf_press(p.I), gf_velx(p.I),
                   gf_vely(p.I), gf_velz(p.I), w_lorentz, Ex, Ey, Ez, Bx, By,
                   Bz);

        /* Write back conserved variables in case they have been adjusted
           or set to NaN                                                    */
        if (rep.adjust_cons)
          cv.scatter(gf_dens(p.I), gf_tau(p.I), Y_e_con, gf_momx(p.I),
                     gf_momy(p.I), gf_momy(p.I), dBx, dBy, dBz);
      });
}
} // namespace GRHydroToyGPU
