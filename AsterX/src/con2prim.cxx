#include <fixmath.hxx>
#include <loop_device.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>
#include "con2prim.h"

namespace AsterX {
using namespace std;
using namespace Loop;

#ifdef __cplusplus
extern "C" {
#endif
/***************************************************************************
2D Newton-Raphson for just-hydro ideal-fluid c2p
------------------------------------
Given an initial guess x[0], x[1] for a root, take ntrial Newton-Raphson steps
until reaching tolf. Called from AsterX_Con2Prim_2DNRNoble. This has been
inspired by NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING (ISBN
0-521-43108-5).
****************************************************************************/
CCTK_DEVICE CCTK_INT mnewt(CCTK_INT ntrial, CCTK_REAL x[], CCTK_REAL Dens,
                           CCTK_REAL tau, CCTK_REAL S_squared,
                           CCTK_REAL gamma) {
  const CCTK_INT n = 2;
  CCTK_REAL tolf = 1E-9;
  CCTK_INT k, i;
  CCTK_REAL fvec[n];
  CCTK_REAL dx[n];
  CCTK_REAL fjac[n][n];

  CCTK_REAL rho, press, dPdx0, dPdx1;
  CCTK_REAL a, b, c, d, detjac_inv;
  CCTK_REAL errf;

  for (k = 1; k <= ntrial; k++) {
    rho = Dens * sqrt(1.0 - x[1]);
    press = (x[0] * (1.0 - x[1]) / rho - 1.0) * rho * (gamma - 1.0) / gamma;
    fvec[0] = x[1] * x[0] * x[0] - S_squared;
    fvec[1] = tau + Dens - x[0] + press;
    dPdx0 = (1.0 - x[1]) * (gamma - 1.0) / gamma;
    dPdx1 = -x[0] * (gamma - 1.0) / gamma;
    fjac[0][0] = 2.0 * x[1] * x[0];
    fjac[0][1] = x[0] * x[0];
    fjac[1][0] = -1.0 + dPdx0;
    fjac[1][1] = dPdx1;
    a = fjac[0][0];
    b = fjac[0][1];
    c = fjac[1][0];
    d = fjac[1][1];
    detjac_inv = 1.0 / (a * d - b * c);
    dx[0] = -detjac_inv * (d * fvec[0] - b * fvec[1]);
    dx[1] = -detjac_inv * (-c * fvec[0] + a * fvec[1]);

    errf = 0.0;
    for (i = 0; i < n; i++) {
      errf += fabs(fvec[i]);
    }
    if (errf <= tolf) {
      return 0;
    }
    for (i = 0; i < n; i++) {
      x[i] += dx[i];
    }
  }
  return 1;
}
#ifdef __cplusplus
}
#endif

/***************************************************************************
2DNRNoble just-hydro ideal-fluid flat-space C2P
------------------------------------
2D-NR Noble scheme for c2p.
Valid for just-hydro, ideal-fluid, and flat-space.
This has been inspired by Siegel+2018.
****************************************************************************/
extern "C" void AsterX_Con2Prim_2DNRNoble(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Con2Prim_2DNRNoble;
  DECLARE_CCTK_PARAMETERS;

  constexpr auto DI = PointDesc::DI;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  constexpr array<int, dim> vertex_centred = {0, 0, 0};
  const GF3D2layout gf_layout_cell(cctkGH, cell_centred);
  const GF3D2layout gf_layout_vertex(cctkGH, vertex_centred);

  const GF3D2<const CCTK_REAL> gf_alp(gf_layout_vertex, alp);
  const GF3D2<const CCTK_REAL> gf_betax(gf_layout_vertex, betax);
  const GF3D2<const CCTK_REAL> gf_betay(gf_layout_vertex, betay);
  const GF3D2<const CCTK_REAL> gf_betaz(gf_layout_vertex, betaz);
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

  // Loop over the entire grid (0 to n-1 points in each direction)
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        CCTK_REAL alp = 0.0;   // lapse
        CCTK_REAL betax = 0.0; // beta^i
        CCTK_REAL betay = 0.0;
        CCTK_REAL betaz = 0.0;
        CCTK_REAL betalx, betaly, betalz; // beta_i

        CCTK_REAL g_up[4][4] = {
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
        };

        CCTK_REAL g_lo[4][4] = {
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
            {0.0, 0.0, 0.0, 0.0},
        };

        for (int dk = 0; dk < 2; ++dk)
          for (int dj = 0; dj < 2; ++dj)
            for (int di = 0; di < 2; ++di) {
              g_lo[1][1] += gf_gxx(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g_lo[1][2] += gf_gxy(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g_lo[1][3] += gf_gxz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g_lo[2][2] += gf_gyy(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g_lo[2][3] += gf_gyz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g_lo[3][3] += gf_gzz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              alp += gf_alp(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              betax += gf_betax(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              betay += gf_betay(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              betaz += gf_betaz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
            }

        alp *= 0.125;
        betax *= 0.125;
        betay *= 0.125;
        betaz *= 0.125;
        g_lo[1][1] *= 0.125;
        g_lo[1][2] *= 0.125;
        g_lo[1][3] *= 0.125;
        g_lo[2][2] *= 0.125;
        g_lo[2][3] *= 0.125;
        g_lo[3][3] *= 0.125;

        betalx = g_lo[1][1] * betax + g_lo[1][2] * betay + g_lo[1][3] * betaz;
        betaly = g_lo[1][2] * betax + g_lo[2][2] * betay + g_lo[2][3] * betaz;
        betalz = g_lo[1][3] * betax + g_lo[2][3] * betay + g_lo[3][3] * betaz;

        g_lo[0][0] =
            -alp * alp + betax * betalx + betay * betaly + betaz * betalz;
        g_lo[0][1] = betalx;
        g_lo[0][2] = betaly;
        g_lo[0][3] = betalz;

        g_lo[1][0] = g_lo[0][1];
        g_lo[2][0] = g_lo[0][2];
        g_lo[3][0] = g_lo[0][3];
        g_lo[2][1] = g_lo[1][2];
        g_lo[3][1] = g_lo[1][3];
        g_lo[3][2] = g_lo[2][3];

        /* Calculate inverse of 4-dim metric */
        // TODO: Move to function?
        CCTK_REAL spatial_detg; // Determinant spatial metric
        CCTK_REAL gamma11, gamma12, gamma13, gamma22, gamma23,
            gamma33; // Inverse components of spatial metric

        spatial_detg = -g_lo[1][3] * g_lo[1][3] * g_lo[2][2] +
                       2 * g_lo[1][2] * g_lo[1][3] * g_lo[2][3] -
                       g_lo[1][1] * g_lo[2][3] * g_lo[2][3] -
                       g_lo[1][2] * g_lo[1][2] * g_lo[3][3] +
                       g_lo[1][1] * g_lo[2][2] * g_lo[3][3];

        gamma11 =
            (-g_lo[2][3] * g_lo[2][3] + g_lo[2][2] * g_lo[3][3]) / spatial_detg;
        gamma12 =
            (g_lo[1][3] * g_lo[2][3] - g_lo[1][2] * g_lo[3][3]) / spatial_detg;
        gamma13 =
            (-g_lo[1][3] * g_lo[2][2] + g_lo[1][2] * g_lo[2][3]) / spatial_detg;
        gamma22 =
            (-g_lo[1][3] * g_lo[1][3] + g_lo[1][1] * g_lo[3][3]) / spatial_detg;
        gamma23 =
            (g_lo[1][2] * g_lo[1][3] - g_lo[1][1] * g_lo[2][3]) / spatial_detg;
        gamma33 =
            (-g_lo[1][2] * g_lo[1][2] + g_lo[1][1] * g_lo[2][2]) / spatial_detg;

        g_up[0][0] = -1.0 / (alp * alp);
        g_up[0][1] = betax / (alp * alp);
        g_up[0][2] = betay / (alp * alp);
        g_up[0][3] = betaz / (alp * alp);
        g_up[1][1] = gamma11 - betax * betax / (alp * alp);
        g_up[1][2] = gamma12 - betax * betay / (alp * alp);
        g_up[1][3] = gamma13 - betax * betaz / (alp * alp);
        g_up[2][2] = gamma22 - betay * betay / (alp * alp);
        g_up[2][3] = gamma23 - betay * betaz / (alp * alp);
        g_up[3][3] = gamma33 - betaz * betaz / (alp * alp);

        g_up[1][0] = g_up[0][1];
        g_up[2][0] = g_up[0][2];
        g_up[3][0] = g_up[0][3];
        g_up[2][1] = g_up[1][2];
        g_up[3][1] = g_up[1][3];
        g_up[3][2] = g_up[2][3];

        CCTK_REAL v1_coord_con, v2_coord_con, v3_coord_con;
        CCTK_REAL v1_coord_cov, v2_coord_cov, v3_coord_cov;
        CCTK_REAL S_squared;
        CCTK_REAL vsq_last, gamma_last;
        CCTK_REAL rho_last, eps_last, Z_last, press_last;
        CCTK_REAL x_2d[2];
        CCTK_INT c2p_failed;
        CCTK_REAL Z, W;

        /* prims from last timestep */
        CCTK_REAL prim[NPRIMS];
        prim[RHO] = gf_rho(p.I);
        // contravariant coordinate velocity:
        v1_coord_con = gf_velx(p.I);
        v2_coord_con = gf_vely(p.I);
        v3_coord_con = gf_velz(p.I);
        // convariant coordinate velocity:
        v1_coord_cov = g_lo[1][1] * v1_coord_con + g_lo[1][2] * v2_coord_con +
                       g_lo[1][3] * v3_coord_con;
        v2_coord_cov = g_lo[2][1] * v1_coord_con + g_lo[2][2] * v2_coord_con +
                       g_lo[2][3] * v3_coord_con;
        v3_coord_cov = g_lo[3][1] * v1_coord_con + g_lo[3][2] * v2_coord_con +
                       g_lo[3][3] * v3_coord_con;
        // covariant Valencia velocity:
        prim[V1_COV] = (v1_coord_cov + g_lo[0][1]) / alp;
        prim[V2_COV] = (v2_coord_cov + g_lo[0][2]) / alp;
        prim[V3_COV] = (v3_coord_cov + g_lo[0][3]) / alp;
        prim[EPS] = gf_eps(p.I);

        /* cons at this timestep */
        CCTK_REAL con[NCONS];
        con[D] = gf_dens(p.I);
        con[S1_cov] = gf_momx(p.I);
        con[S2_cov] = gf_momy(p.I);
        con[S3_cov] = gf_momz(p.I);
        con[TAU] = gf_tau(p.I);

        /* calculate W from last timestep and use for guess */
        vsq_last = prim[V1_COV] *
                   (g_up[1][1] * prim[V1_COV] + g_up[1][2] * prim[V2_COV] +
                    g_up[1][3] * prim[V3_COV]);
        vsq_last += prim[V2_COV] *
                    (g_up[2][1] * prim[V1_COV] + g_up[2][2] * prim[V2_COV] +
                     g_up[2][3] * prim[V3_COV]);
        vsq_last += prim[V3_COV] *
                    (g_up[3][1] * prim[V1_COV] + g_up[3][2] * prim[V2_COV] +
                     g_up[3][3] * prim[V3_COV]);
        if ((vsq_last < 0.) && (fabs(vsq_last) < 1.0e-13)) {
          vsq_last = fabs(vsq_last);
        }
        gamma_last = 1. / sqrt(1. - vsq_last);

        /* calculate S_squared */
        S_squared =
            con[S1_cov] * (g_up[1][1] * con[S1_cov] + g_up[1][2] * con[S2_cov] +
                           g_up[1][3] * con[S3_cov]);
        S_squared +=
            con[S2_cov] * (g_up[1][2] * con[S1_cov] + g_up[2][2] * con[S2_cov] +
                           g_up[2][3] * con[S3_cov]);
        S_squared +=
            con[S3_cov] * (g_up[1][3] * con[S1_cov] + g_up[2][3] * con[S2_cov] +
                           g_up[3][3] * con[S3_cov]);
        if ((S_squared < 0.) && (fabs(S_squared) < 1.0e-13)) {
          S_squared = fabs(S_squared);
        }

        /* calculate rho from D and gamma */
        rho_last =
            con[D] / gamma_last; // rho consistent with con[D] should be better
                                 // guess than rho from last timestep
        eps_last = prim[EPS];
        press_last = (gamma - 1.0) * eps_last * rho_last;
        Z_last = (rho_last + eps_last * rho_last + press_last) * gamma_last *
                 gamma_last; // Z = rho h W^2

        /* initialize unknowns for c2p, Z and vsq: */
        x_2d[0] = fabs(Z_last);
        x_2d[1] = vsq_last;

        /* start recovery */
        c2p_failed = mnewt(max_iter, x_2d, con[D], con[TAU], S_squared, gamma);

        /* Calculate primitives from Z and W */
        Z = x_2d[0];
        W = 1.0 / sqrt(1.0 - x_2d[1]);
        if (c2p_failed || (Z != Z) || (W != W)) {
          // If c2p fails, reset prims
          // TODO: set counter and error message
          prim[RHO] = rho_last;
          prim[V1_COV] = 0.0;
          prim[V2_COV] = 0.0;
          prim[V3_COV] = 0.0;
          prim[EPS] = eps_last;
          assert(0); // Terminate?
        } else {
          prim[RHO] = con[D] / W;
          prim[V1_COV] = con[S1_cov] / Z;
          prim[V2_COV] = con[S2_cov] / Z;
          prim[V3_COV] = con[S3_cov] / Z;
          prim[EPS] = (Z * (1.0 - x_2d[1]) / prim[RHO] - 1.0) / gamma;
        }
        gf_rho(p.I) = prim[0];
        v1_coord_cov = alp * prim[V1_COV] - g_lo[0][1];
        v2_coord_cov = alp * prim[V2_COV] - g_lo[0][2];
        v3_coord_cov = alp * prim[V3_COV] - g_lo[0][3];
        gf_velx(p.I) = g_up[1][1] * v1_coord_cov + g_up[1][2] * v2_coord_cov +
                       g_up[1][3] * v3_coord_cov;
        gf_vely(p.I) = g_up[1][2] * v1_coord_cov + g_up[2][2] * v2_coord_cov +
                       g_up[1][3] * v3_coord_cov;
        gf_velz(p.I) = g_up[1][3] * v1_coord_cov + g_up[1][3] * v2_coord_cov +
                       g_up[3][3] * v3_coord_cov;
        gf_eps(p.I) = prim[EPS];
        gf_press(p.I) = (gamma - 1.0) * prim[4] * prim[0];
      });
}

} // namespace AsterX
