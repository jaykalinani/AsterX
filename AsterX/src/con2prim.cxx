#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include <cmath>

#include "con2primIdealFluid.hxx"
namespace AsterX {
using namespace std;
using namespace Loop;

/***************************************************************************
2DNRNoble C2P
------------------------------------
2D-NR Noble scheme for c2p.
Sources: Noble+2006, Section 3.1 of Siegel+2018,
NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
****************************************************************************/
template <typename typeEoS>
CCTK_HOST CCTK_DEVICE void Con2Prim_2DNRNoble(
    CCTK_INT max_iter, CCTK_REAL tolf,
    typeEoS &plasma)
{ // Send con2primFactory object as reference to modify it,
  // and because we can not instantiate abstract class

    /* get Lorentz factor seed, calculated by constructor */
    CCTK_REAL W = plasma.WLorentz_Seed;

    /* get pressure seed */
    plasma.get_Press_Seed();

    /* get Z seed */
    plasma.get_Z_Seed();

    /* initialize unknowns for c2p, Z and vsq: */
    CCTK_REAL x[2];
    CCTK_REAL x_old[2];
    x[0] = fabs(plasma.Z_Seed);
    x[1] = (-1.0 + W * W) / (W * W);

    /* initialize old values */
    x_old[0] = x[0];
    x_old[1] = x[1];

    /* Start Recovery with 2D NR Solver */
    const CCTK_INT n = 2;
    const CCTK_REAL dv = (1. - 1.e-15);
    CCTK_REAL fvec[n];
    CCTK_REAL dx[n];
    CCTK_REAL fjac[n][n];

    CCTK_REAL detjac_inv;
    CCTK_REAL errf;
    plasma.Failed_2DNRNoble = 1;
    CCTK_INT k;
    for (k = 1; k <= max_iter; k++)
    {
        fvec[0] = plasma.get_2DNRNoble_f0(x[0], x[1]);
        fvec[1] = plasma.get_2DNRNoble_f1(x[0], x[1]);
        fjac[0][0] = plasma.get_2DNRNoble_df0dZ(x[0], x[1]);
        fjac[0][1] = plasma.get_2DNRNoble_df0dVsq(x[0], x[1]);
        fjac[1][0] = plasma.get_2DNRNoble_df1dZ(x[0], x[1]);
        fjac[1][1] = plasma.get_2DNRNoble_df1dVsq(x[0], x[1]);
        detjac_inv = 1.0 / (fjac[0][0] * fjac[1][1] - fjac[0][1] * fjac[1][0]);
        dx[0] = -detjac_inv * (fjac[1][1] * fvec[0] - fjac[0][1] * fvec[1]);
        dx[1] = -detjac_inv * (-fjac[1][0] * fvec[0] + fjac[0][0] * fvec[1]);

        errf = 0.0;
        for (CCTK_INT i = 0; i < n; i++)
        {
            errf += fabs(fvec[i]);
        }
        if (errf <= tolf)
        {
            plasma.Failed_2DNRNoble = 0;
            break;
        }

        /* save old values before calculating the new */
        x_old[0] = x[0];
        x_old[1] = x[1];

        for (CCTK_INT i = 0; i < n; i++)
        {
            x[i] += dx[i];
        }
    }
    plasma.Nit_2DNRNoble = k;

  /* make sure that the new x[] is physical */
  if (x[0] < 0.0) {
    x[0] = fabs(x[0]);
  } else {
    if (x[0] > 1e20) {
      x[0] = x_old[0];
    }
  }

  if (x[1] < 0.0) {
    x[1] = 0.0;
  } else {
    if (x[1] >= 1.0) {
      x[1] = dv;
    }
  }

  /* Calculate primitives from Z and W */
  plasma.Z_Sol = x[0];
  plasma.vsq_Sol = x[1];
}

///

template <typename typeEoS> void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the interior of the grid
  cctk_grid.loop_int_device<
      1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    CCTK_REAL g_up[4][4];
    CCTK_REAL g_lo[4][4];

    /* Get covariant metric */
    g_lo[1][1] = calc_avg_v2c(gxx, p);
    g_lo[1][2] = calc_avg_v2c(gxy, p);
    g_lo[1][3] = calc_avg_v2c(gxz, p);
    g_lo[2][2] = calc_avg_v2c(gyy, p);
    g_lo[2][3] = calc_avg_v2c(gyz, p);
    g_lo[3][3] = calc_avg_v2c(gzz, p);
    CCTK_REAL lapse = calc_avg_v2c(alp, p);
    CCTK_REAL betax_up = calc_avg_v2c(betax, p);
    CCTK_REAL betay_up = calc_avg_v2c(betay, p);
    CCTK_REAL betaz_up = calc_avg_v2c(betaz, p);

    // beta_lo
    g_lo[0][1] =
        g_lo[1][1] * betax_up + g_lo[1][2] * betay_up + g_lo[1][3] * betaz_up;
    g_lo[0][2] =
        g_lo[1][2] * betax_up + g_lo[2][2] * betay_up + g_lo[2][3] * betaz_up;
    g_lo[0][3] =
        g_lo[1][3] * betax_up + g_lo[2][3] * betay_up + g_lo[3][3] * betaz_up;

    g_lo[0][0] = -lapse * lapse + betax_up * g_lo[0][1] +
                 betay_up * g_lo[0][2] + betaz_up * g_lo[0][3];

    // symmetric components
    g_lo[1][0] = g_lo[0][1];
    g_lo[2][0] = g_lo[0][2];
    g_lo[3][0] = g_lo[0][3];
    g_lo[2][1] = g_lo[1][2];
    g_lo[3][1] = g_lo[1][3];
    g_lo[3][2] = g_lo[2][3];

    /* Calculate inverse of 4-dim metric */
    CCTK_REAL gamma11, gamma12, gamma13, gamma22, gamma23,
        gamma33; // Inverse components of spatial metric

    CCTK_REAL spatial_detg = -g_lo[1][3] * g_lo[1][3] * g_lo[2][2] +
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

    g_up[0][0] = -1.0 / (lapse * lapse);
    g_up[0][1] = betax_up / (lapse * lapse);
    g_up[0][2] = betay_up / (lapse * lapse);
    g_up[0][3] = betaz_up / (lapse * lapse);
    g_up[1][1] = gamma11 - betax_up * betax_up / (lapse * lapse);
    g_up[1][2] = gamma12 - betax_up * betay_up / (lapse * lapse);
    g_up[1][3] = gamma13 - betax_up * betaz_up / (lapse * lapse);
    g_up[2][2] = gamma22 - betay_up * betay_up / (lapse * lapse);
    g_up[2][3] = gamma23 - betay_up * betaz_up / (lapse * lapse);
    g_up[3][3] = gamma33 - betaz_up * betaz_up / (lapse * lapse);

    g_up[1][0] = g_up[0][1];
    g_up[2][0] = g_up[0][2];
    g_up[3][0] = g_up[0][3];
    g_up[2][1] = g_up[1][2];
    g_up[3][1] = g_up[1][3];
    g_up[3][2] = g_up[2][3];

    // Set to atmosphere if dens below the threshold
    if (dens(p.I) <= sqrt(spatial_detg) * (rho_abs_min * (1.0 + atmo_tol))) {
      saved_rho(p.I) = rho_abs_min;
      saved_velx(p.I) = 0.0;
      saved_vely(p.I) = 0.0;
      saved_velz(p.I) = 0.0;
      press(p.I) = poly_K * pow(saved_rho(p.I), gamma);
      saved_eps(p.I) = press(p.I) / ((gamma - 1.0) * saved_rho(p.I));

      dens(p.I) = sqrt(spatial_detg) * saved_rho(p.I);
      momx(p.I) = 0.0;
      momy(p.I) = 0.0;
      momz(p.I) = 0.0;
      // need to compute bs2; setting here to 0.0
      CCTK_REAL bs2 = 0.0;
      tau(p.I) = sqrt(spatial_detg) *
                     (saved_rho(p.I) * (1 + saved_eps(p.I)) + 0.5 * bs2) -
                 dens(p.I);
    }

    CCTK_REAL cons[NCONS];
    cons[D] = dens(p.I) / sqrt(spatial_detg);
    cons[S1_COV] = momx(p.I) / sqrt(spatial_detg);
    cons[S2_COV] = momy(p.I) / sqrt(spatial_detg);
    cons[S3_COV] = momz(p.I) / sqrt(spatial_detg);
    cons[TAU] = tau(p.I) / sqrt(spatial_detg);
    cons[B1] = dBx(p.I) / sqrt(spatial_detg);
    cons[B2] = dBy(p.I) / sqrt(spatial_detg);
    cons[B3] = dBz(p.I) / sqrt(spatial_detg);

    CCTK_REAL prims[NPRIMS];
    prims[RHO] = saved_rho(p.I);
    prims[V1_CON] = saved_velx(p.I);
    prims[V2_CON] = saved_vely(p.I);
    prims[V3_CON] = saved_velz(p.I);
    prims[EPS] = saved_eps(p.I);
    prims[B1] = cons[B1];
    prims[B2] = cons[B2]; 
    prims[B3] = cons[B3]; 


    // Construct con2primFactory object:
    typeEoS plasma_0(gamma, cons, prims, g_lo, g_up);

    // 1) Try 2DNRNoble
    Con2Prim_2DNRNoble(max_iter, c2p_tol, plasma_0);

    CCTK_INT set_to_atmo = 0;

    if (plasma_0.Failed_2DNRNoble) {
      if (debug_mode) {
        printf(
            "WARNING: "
            "2DNRNoble failed. Printing cons and saved prims before set to "
            "atmo: \n"
            "cctk_iteration = %i \n "
            "x, y, z = %26.16e, %26.16e, %26.16e \n "
            "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
            "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
            "dBy = %26.16e \n dBz = %26.16e \n "
            "saved_rho = %26.16e \n saved_eps = %26.16e \n press= %26.16e \n "
            "saved_velx = %26.16e \n saved_vely = %26.16e \n saved_velz = "
            "%26.16e \n "
            "Bvecx = %26.16e \n Bvecy = %26.16e \n "
            "Bvecz = %26.16e \n "
            "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
            cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I), momx(p.I),
            momy(p.I), momz(p.I), dBx(p.I), dBy(p.I), dBz(p.I), saved_rho(p.I),
            saved_eps(p.I), press(p.I), saved_velx(p.I), saved_vely(p.I),
            saved_velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I), Avec_x(p.I),
            Avec_y(p.I), Avec_z(p.I));
      }

      set_to_atmo = 1;
      // CCTK_VWARN(CCTK_WARN_ALERT, "C2P failed, triggering assert(0)..");
      // assert(0);
    } else {
      plasma_0.WZ2Prim();
      rho(p.I) = plasma_0.PrimitiveVars[RHO];
      velx(p.I) = plasma_0.PrimitiveVars[V1_CON];
      vely(p.I) = plasma_0.PrimitiveVars[V2_CON];
      velz(p.I) = plasma_0.PrimitiveVars[V3_CON];
      eps(p.I) = plasma_0.PrimitiveVars[EPS];
      press(p.I) = (gamma - 1.0) * eps(p.I) * rho(p.I);

      if (rho(p.I) <= rho_abs_min * (1 + atmo_tol)) {
        set_to_atmo = 1;
      }
    }

    if (set_to_atmo) {
      rho(p.I) = rho_abs_min;
      velx(p.I) = 0.0;
      vely(p.I) = 0.0;
      velz(p.I) = 0.0;
      press(p.I) = poly_K * pow(rho(p.I), gamma);
      eps(p.I) = press(p.I) / ((gamma - 1.0) * rho(p.I));

      dens(p.I) = sqrt(spatial_detg) * rho(p.I);
      momx(p.I) = 0.0;
      momy(p.I) = 0.0;
      momz(p.I) = 0.0;
      // need to compute bs2; setting here to 0.0
      CCTK_REAL bs2 = 0.0;
      tau(p.I) = sqrt(spatial_detg) * (rho(p.I) * (1 + eps(p.I)) + 0.5 * bs2) -
                 dens(p.I);

      if (debug_mode) {
        printf("WARNING: "
               "Printing cons and prims after set to atmo when dens below "
               "threshold: \n"
               "cctk_iteration = %i \n "
               "x, y, z = %26.16e, %26.16e, %26.16e \n "
               "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
               "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
               "dBy = %26.16e \n dBz = %26.16e \n "
               "rho = %26.16e \n eps = %26.16e \n press= %26.16e \n "
               "velx = %26.16e \n vely = %26.16e \n velz = %26.16e \n "
               "Bvecx = %26.16e \n Bvecy = %26.16e \n "
               "Bvecz = %26.16e \n "
               "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
               cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I), momx(p.I),
               momy(p.I), momz(p.I), dBx(p.I), dBy(p.I), dBz(p.I), rho(p.I),
               eps(p.I), press(p.I), velx(p.I), vely(p.I), velz(p.I),
               Bvecx(p.I), Bvecy(p.I), Bvecz(p.I), Avec_x(p.I), Avec_y(p.I),
               Avec_z(p.I));
      }
    }

    Bvecx(p.I) = cons[B1];
    Bvecy(p.I) = cons[B2];
    Bvecz(p.I) = cons[B3];
    saved_rho(p.I) = rho(p.I);
    saved_velx(p.I) = velx(p.I);
    saved_vely(p.I) = vely(p.I);
    saved_velz(p.I) = velz(p.I);
    saved_eps(p.I) = eps(p.I);
  }); // Loop
} // AsterX_Con2Prim_2DNRNoble

/***************************************************************************
 * AsterX_Con2Prim
 * ---
 *  Routines implemented:
 *   1) 2DNRNoble
 *
 *   Based on con2primFactory (https://github.com/fedelopezar/con2primFactory)
 *   ****************************************************************************/
extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS) {
  if (1) { // Use this if for idealFluid/tabeos
    AsterX_Con2Prim_typeEoS<idealFluid>(CCTK_PASS_CTOC);
    // CCTK_PASS_CTOC == cctkGH, and more. Preferred over just cctkGH.
  }
}

} // namespace AsterX
