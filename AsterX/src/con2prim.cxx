#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>

#include <cmath>

#include "con2prim.hxx"
namespace AsterX
{
  using namespace std;
  using namespace Loop;

  /* Constructor */
  CCTK_HOST CCTK_DEVICE idealFluid::idealFluid(CCTK_REAL gamma, CCTK_REAL (&cons)[NCONS], CCTK_REAL (&prim)[NPRIMS], CCTK_REAL (&g_lo)[4][4], CCTK_REAL (&g_up)[4][4])
  {
    GammaIdealFluid = gamma;

    PrimitiveVarsSeed[RHO] = prim[RHO];
    PrimitiveVarsSeed[V1_CON] = prim[V1_CON];
    PrimitiveVarsSeed[V2_CON] = prim[V2_CON];
    PrimitiveVarsSeed[V3_CON] = prim[V3_CON];
    PrimitiveVarsSeed[EPS] = prim[EPS];
    PrimitiveVarsSeed[B1] = prim[B1];
    PrimitiveVarsSeed[B2] = prim[B2];
    PrimitiveVarsSeed[B3] = prim[B3];

    ConservedVars[D] = cons[D];
    ConservedVars[S1_COV] = cons[S1_COV];
    ConservedVars[S2_COV] = cons[S2_COV];
    ConservedVars[S3_COV] = cons[S3_COV];
    ConservedVars[EPS] = cons[EPS];
    ConservedVars[B1] = cons[B1];
    ConservedVars[B2] = cons[B2];
    ConservedVars[B3] = cons[B3];

    gcov[TT] = g_lo[0][0];
    gcov[TX] = g_lo[0][1];
    gcov[TY] = g_lo[0][2];
    gcov[TZ] = g_lo[0][3];
    gcov[XX] = g_lo[1][1];
    gcov[XY] = g_lo[1][2];
    gcov[XZ] = g_lo[1][3];
    gcov[YY] = g_lo[2][2];
    gcov[YZ] = g_lo[2][3];
    gcov[ZZ] = g_lo[3][3];

    gcon[TT] = g_up[0][0];
    gcon[TX] = g_up[0][1];
    gcon[TY] = g_up[0][2];
    gcon[TZ] = g_up[0][3];
    gcon[XX] = g_up[1][1];
    gcon[XY] = g_up[1][2];
    gcon[XZ] = g_up[1][3];
    gcon[YY] = g_up[2][2];
    gcon[YZ] = g_up[2][3];
    gcon[ZZ] = g_up[3][3];

    CCTK_REAL alp = sqrt(-1. / gcon[TT]);

    /* B^i S_i */
    BiSi = prim[B1] * ConservedVars[S1_COV] + prim[B2] * ConservedVars[S2_COV] + prim[B3] * ConservedVars[S3_COV];

    /* Seed Lorentz factor */
    // covariant Valencia velocity:
    CCTK_REAL v1_cov = gcov[XX] * PrimitiveVarsSeed[V1_CON] + gcov[XY] * PrimitiveVarsSeed[V2_CON] +
                       gcov[XZ] * PrimitiveVarsSeed[V3_CON];
    CCTK_REAL v2_cov = gcov[XY] * PrimitiveVarsSeed[V1_CON] + gcov[YY] * PrimitiveVarsSeed[V2_CON] +
                       gcov[YZ] * PrimitiveVarsSeed[V3_CON];
    CCTK_REAL v3_cov = gcov[XZ] * PrimitiveVarsSeed[V1_CON] + gcov[YZ] * PrimitiveVarsSeed[V2_CON] +
                       gcov[ZZ] * PrimitiveVarsSeed[V3_CON];

    CCTK_REAL vsq = v1_cov * PrimitiveVarsSeed[V1_CON] + v2_cov * PrimitiveVarsSeed[V2_CON] + v3_cov * PrimitiveVarsSeed[V3_CON];

    if ((vsq < 0.) && (fabs(vsq) < 1.0e-13))
    {
      vsq = fabs(vsq);
    }

    W_Seed = 1. / sqrt(1. - vsq);

    // Bsq and bsq:
    CCTK_REAL B1_cov = gcov[XX] * prim[B1] + gcov[XY] * prim[B2] + gcov[XZ] * prim[B3];
    CCTK_REAL B2_cov = gcov[XY] * prim[B1] + gcov[YY] * prim[B2] + gcov[YZ] * prim[B3];
    CCTK_REAL B3_cov = gcov[XZ] * prim[B1] + gcov[YZ] * prim[B2] + gcov[ZZ] * prim[B3];

    Bsq = B1_cov * prim[B1] + B2_cov * prim[B2] + B3_cov * prim[B3];

    CCTK_REAL bt = W_Seed / alp * (prim[B1] * v1_cov + prim[B2] * v2_cov + prim[B3] * v3_cov);

    bsq = (Bsq + (alp * bt) * (alp * bt)) / (W_Seed * W_Seed);

    /* update rho seed from D and gamma */
    // rho consistent with con[D] should be better guess than rho from last timestep
    PrimitiveVarsSeed[0] = ConservedVars[D] / W_Seed;
  }

  /* Called by 2dNRNoble */
  CCTK_HOST CCTK_DEVICE void idealFluid::get_Ssq_Exact()
  {

    /* calculate S_squared */
    Ssq =
        ConservedVars[S1_COV] * (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] +
                                 gcon[XZ] * ConservedVars[S3_COV]);
    Ssq +=
        ConservedVars[S2_COV] * (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] +
                                 gcon[YZ] * ConservedVars[S3_COV]);
    Ssq +=
        ConservedVars[S3_COV] * (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] +
                                 gcon[ZZ] * ConservedVars[S3_COV]);
    if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13))
    {
      Ssq = fabs(Ssq);
    }
  }

  CCTK_HOST CCTK_DEVICE void idealFluid::get_Press_Seed()
  {
    Press_Seed = PrimitiveVarsSeed[RHO] * PrimitiveVarsSeed[EPS] * (GammaIdealFluid - 1.0);
  }

  CCTK_HOST CCTK_DEVICE void idealFluid::get_Z_Seed()
  {
    Z_Seed = (PrimitiveVarsSeed[RHO] + PrimitiveVarsSeed[EPS] * PrimitiveVarsSeed[RHO] + Press_Seed) * W_Seed * W_Seed;
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    return (Vsq * (Bsq + Z) * (Bsq + Z) - (BiSi * BiSi * (Bsq + 2.0 * Z)) / (Z * Z) - Ssq);
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    CCTK_REAL Press = get_Press_funcZVsq(Z, Vsq);
    return ConservedVars[TAU] + ConservedVars[D] - Bsq / 2.0 * (1 + Vsq) + BiSi * BiSi / 2.0 / (Z * Z) - Z + Press;
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    return (2.0 * Vsq * (Bsq + Z) - 2.0 * BiSi * BiSi / (Z * Z) + 2.0 * BiSi * BiSi * (Bsq + 2.0 * Z) / (Z * Z * Z));
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    return (Bsq + Z) * (Bsq + Z);
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    return (-BiSi * BiSi / (Z * Z * Z) - 1.0 + get_dPdZ_funcZVsq(Z, Vsq));
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_2DNRNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    return (-Bsq / 2.0 + get_dPdVsq_funcZVsq(Z, Vsq));
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    return ((Z * (1.0 - Vsq) - ConservedVars[D] * sqrt(1.0 - Vsq)) * (GammaIdealFluid - 1.0) / (GammaIdealFluid));
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    return ((1.0 - Vsq) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
  }

  CCTK_HOST CCTK_DEVICE CCTK_REAL idealFluid::get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq)
  {
    return ((-Z + ConservedVars[D] / (2.0 * sqrt(1.0 - Vsq))) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
  }

  CCTK_HOST CCTK_DEVICE void idealFluid::WZ2Prim()
  {
    CCTK_REAL W_Sol = 1.0 / sqrt(1.0 - vsq_Sol);
    PrimitiveVars[RHO] = ConservedVars[D] / W_Sol;
    CCTK_REAL alp = sqrt(-1. / gcon[TT]);

    PrimitiveVars[V1_CON] = (gcon[XX] * ConservedVars[S1_COV] + gcon[XY] * ConservedVars[S2_COV] + gcon[XZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    PrimitiveVars[V1_CON] += BiSi * PrimitiveVars[B1] / (Z_Sol * (Z_Sol + Bsq));

    PrimitiveVars[V2_CON] = (gcon[XY] * ConservedVars[S1_COV] + gcon[YY] * ConservedVars[S2_COV] + gcon[YZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    PrimitiveVars[V2_CON] += BiSi * PrimitiveVars[B2] / (Z_Sol * (Z_Sol + Bsq));

    PrimitiveVars[V3_CON] = (gcon[XZ] * ConservedVars[S1_COV] + gcon[YZ] * ConservedVars[S2_COV] + gcon[ZZ] * ConservedVars[S3_COV]) / (Z_Sol + Bsq);
    PrimitiveVars[V3_CON] += BiSi * PrimitiveVars[B3] / (Z_Sol * (Z_Sol + Bsq));

    PrimitiveVars[EPS] = (Z_Sol * (1. - vsq_Sol) / PrimitiveVars[RHO] - 1.0) / GammaIdealFluid;
    PrimitiveVars[B1] = ConservedVars[B1];
    PrimitiveVars[B2] = ConservedVars[B2];
    PrimitiveVars[B3] = ConservedVars[B3];
  }

  /* Destructor */
  CCTK_HOST CCTK_DEVICE idealFluid::~idealFluid()
  {
    // How to destruct properly a vector?
  }

  /***************************************************************************
2DNRNoble C2P
------------------------------------
2D-NR Noble scheme for c2p.
Sources: Noble+2006, Section 3.1 of Siegel+2018, 
NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
****************************************************************************/
  template <typename typeEoS>
  CCTK_HOST CCTK_DEVICE void Con2Prim_2DNRNoble(CCTK_INT max_iter, CCTK_REAL tolf, typeEoS &plasma)
  { // Send con2primFactory object as reference to modify it, and because we can not instantiate abstract class

    /* get Lorentz factor seed, calculated by constructor */
    CCTK_REAL W = plasma.W_Seed;

    /* get Ssq from cons (exact) */
    plasma.get_Ssq_Exact();

    /* get pressure seed */
    plasma.get_Press_Seed();

    /* get Z seed */
    plasma.get_Z_Seed();

    /* initialize unknowns for c2p, Z and vsq: */
    CCTK_REAL x[2];
    x[0] = fabs(plasma.Z_Seed);
    x[1] = (-1.0 + W * W) / (W * W);

    /* Start Recovery with 2D NR Solver */
    const CCTK_INT n = 2;
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
      for (CCTK_INT i = 0; i < n; i++)
      {
        x[i] += dx[i];
      }
    }
    plasma.Nit_2DNRNoble = k;

    /* Calculate primitives from Z and W */
    plasma.Z_Sol = x[0];
    plasma.vsq_Sol = x[1];
  }

  ///

  template <typename typeEoS>
  void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS)
  {
    DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
    DECLARE_CCTK_PARAMETERS;

    // Loop over the interior of the grid
    cctk_grid.loop_int_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE
        {
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
          g_lo[0][1] = g_lo[1][1] * betax_up + g_lo[1][2] * betay_up +
                       g_lo[1][3] * betaz_up;
          g_lo[0][2] = g_lo[1][2] * betax_up + g_lo[2][2] * betay_up +
                       g_lo[2][3] * betaz_up;
          g_lo[0][3] = g_lo[1][3] * betax_up + g_lo[2][3] * betay_up +
                       g_lo[3][3] * betaz_up;

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
          g_up[0][1] = g_lo[0][1] / (lapse * lapse);
          g_up[0][2] = g_lo[0][2] / (lapse * lapse);
          g_up[0][3] = g_lo[0][3] / (lapse * lapse);
          g_up[1][1] = gamma11 - g_lo[0][1] * g_lo[0][1] / (lapse * lapse);
          g_up[1][2] = gamma12 - g_lo[0][1] * g_lo[0][2] / (lapse * lapse);
          g_up[1][3] = gamma13 - g_lo[0][1] * g_lo[0][3] / (lapse * lapse);
          g_up[2][2] = gamma22 - g_lo[0][2] * g_lo[0][2] / (lapse * lapse);
          g_up[2][3] = gamma23 - g_lo[0][2] * g_lo[0][3] / (lapse * lapse);
          g_up[3][3] = gamma33 - g_lo[0][3] * g_lo[0][3] / (lapse * lapse);

          g_up[1][0] = g_up[0][1];
          g_up[2][0] = g_up[0][2];
          g_up[3][0] = g_up[0][3];
          g_up[2][1] = g_up[1][2];
          g_up[3][1] = g_up[1][3];
          g_up[3][2] = g_up[2][3];

          CCTK_REAL cons[NCONS];
          cons[D] = dens(p.I);
          cons[S1_COV] = momx(p.I);
          cons[S2_COV] = momy(p.I);
          cons[S3_COV] = momz(p.I);
          cons[TAU] = tau(p.I);
          cons[B1] = Bvecx(p.I);
          cons[B2] = Bvecy(p.I);
          cons[B3] = Bvecz(p.I);

          CCTK_REAL prims[NPRIMS];
          prims[RHO] = saved_rho(p.I);
          prims[V1_CON] = saved_velx(p.I);
          prims[V2_CON] = saved_vely(p.I);
          prims[V3_CON] = saved_velz(p.I);
          prims[EPS] = saved_eps(p.I);
          prims[B1] = saved_Bvecx(p.I);
          prims[B2] = saved_Bvecy(p.I);
          prims[B3] = saved_Bvecz(p.I);

          // Construct con2primFactory object:
          typeEoS plasma_0(gamma, cons, prims, g_lo, g_up);
          // 1) Try 2DNRNoble
          Con2Prim_2DNRNoble( max_iter, c2p_tol, plasma_0);

          if (plasma_0.Failed_2DNRNoble)
          {
            // If c2p fails, reset prims
            rho(p.I) = prims[RHO];
            velx(p.I) = 0.0;
            vely(p.I) = 0.0;
            velz(p.I) = 0.0;
            eps(p.I) = prims[EPS];
            assert(0); // Terminate?
          }
          else
          {
            plasma_0.WZ2Prim();
            rho(p.I) = plasma_0.PrimitiveVars[RHO];
            velx(p.I) = plasma_0.PrimitiveVars[V1_CON];
            vely(p.I) = plasma_0.PrimitiveVars[V2_CON];
            velz(p.I) = plasma_0.PrimitiveVars[V3_CON];
            eps(p.I) = plasma_0.PrimitiveVars[EPS];
            press(p.I) = (gamma - 1.0) * eps(p.I) * rho(p.I);
          }

          saved_rho(p.I) = rho(p.I);
          saved_velx(p.I) = velx(p.I);
          saved_vely(p.I) = vely(p.I);
          saved_velz(p.I) = velz(p.I);
          saved_eps(p.I) = eps(p.I);
          saved_Bvecx(p.I) = Bvecx(p.I);
          saved_Bvecy(p.I) = Bvecy(p.I);
          saved_Bvecz(p.I) = Bvecz(p.I);
        }); // Loop
  }         // AsterX_Con2Prim_2DNRNoble

  /***************************************************************************
 * AsterX_Con2Prim
 * ---
 *  Routines implemented:
 *   1) 2DNRNoble
 *
 *   Based on con2primFactory (https://github.com/fedelopezar/con2primFactory)
 *   ****************************************************************************/
  extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS)
  {
    if (1)
    { // Use this if for idealFluid/tabeos
      AsterX_Con2Prim_typeEoS<idealFluid>(CCTK_PASS_CTOC);
    }
  }

} // namespace AsterX
