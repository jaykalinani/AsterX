#ifndef C2P_2DNOBLE_HXX
#define C2P_2DNOBLE_HXX

#include "c2p.hxx"

namespace Con2PrimFactory {

class c2p_2DNoble : public c2p {
public:
  /* Some attributes */
  CCTK_REAL GammaIdealFluid;
  CCTK_INT Failed_2DNoble;
  CCTK_INT Nit_2DNoble;

  /* Constructor */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble(
      EOSType &eos_th, CCTK_INT maxIter, CCTK_REAL tol);
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Ssq_Exact(vec<CCTK_REAL, 3> &mom, const smat<CCTK_REAL, 3> &gup) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up, const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec, vec<CCTK_REAL, 3> &mom) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 2>
  get_WLorentz_bsq_Seeds(vec<CCTK_REAL, 3> &B_up, vec<CCTK_REAL, 3> &v_up,
                         const smat<CCTK_REAL, 3> &glo) const;

  /* Called by 2DNoble */
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Z_Seed(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL press,
             CCTK_REAL w_lor) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq, cons_vars &cv) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq, cons_vars &cv) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_2DNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL Ssq, CCTK_REAL Bsq,
                 CCTK_REAL BiSi) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_2DNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL Bsq, CCTK_REAL BiSi,
                 cons_vars &cv) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_2DNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL Bsq,
                    CCTK_REAL BiSi) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_2DNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Bsq) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_2DNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL BiSi) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_2DNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL Bsq,
                      cons_vars &cv) const;
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  WZ2Prim(CCTK_REAL Z_Sol, CCTK_REAL vsq_Sol, CCTK_REAL Bsq, CCTK_REAL BiSi,
          EOSType &eos_th, prim_vars &pv, cons_vars &cv,
          const smat<CCTK_REAL, 3> &gup) const;
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds, cons_vars &cv,
        const smat<CCTK_REAL, 3> &glo, CCTK_INT &c2p_succeeded) const;

  /* Destructor */
  CCTK_HOST CCTK_DEVICE ~c2p_2DNoble();
};

/* Constructor */
template <typename EOSType>
CCTK_HOST
    CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble::c2p_2DNoble(
        EOSType &eos_th, CCTK_INT maxIter, CCTK_REAL tol) {
  GammaIdealFluid = eos_th.gamma;
  maxIterations = maxIter;
  tolerance = tol;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Ssq_Exact(vec<CCTK_REAL, 3> &mom_low,
                           const smat<CCTK_REAL, 3> &gup) const {
  vec<CCTK_REAL, 3> mom_up = calc_contraction(gup, mom_low);
  CCTK_REAL Ssq = calc_contraction(mom_low, mom_up);

  if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13)) {
    Ssq = fabs(Ssq);
  }
  return Ssq;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up,
                           const smat<CCTK_REAL, 3> &glo) const {
  vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);
  return calc_contraction(B_low, B_up); // Bsq
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec,
                            vec<CCTK_REAL, 3> &mom) const {
  return calc_contraction(mom, Bvec); // BiS^i
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 2>
c2p_2DNoble::get_WLorentz_bsq_Seeds(vec<CCTK_REAL, 3> &B_up,
                                    vec<CCTK_REAL, 3> &v_up,
                                    const smat<CCTK_REAL, 3> &glo) const {
  vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);
  CCTK_REAL vsq = calc_contraction(v_low, v_up);
  CCTK_REAL VdotB = calc_contraction(v_low, B_up);
  CCTK_REAL VdotBsq = VdotB * VdotB;
  CCTK_REAL Bsq = get_Bsq_Exact(B_up, glo);

  if ((vsq < 0.) && (fabs(vsq) < 1.0e-13)) {
    vsq = fabs(vsq);
  }

  CCTK_REAL w_lor = 1. / sqrt(1. - vsq);
  CCTK_REAL bsq = ((Bsq) / (w_lor * w_lor)) + VdotBsq;
  vec<CCTK_REAL, 2> w_bsq{w_lor, bsq};

  return w_bsq; //{w_lor, bsq}
}

/* Called by 2dNRNoble */

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Z_Seed(CCTK_REAL rho, CCTK_REAL eps, CCTK_REAL press,
                        CCTK_REAL w_lor) const {
  return (rho + eps * rho + press) * w_lor * w_lor;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_Press_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq,
                                cons_vars &cv) const {
  return ((Z * (1.0 - Vsq) - cv.dens * sqrt(1.0 - Vsq)) *
          (GammaIdealFluid - 1.0) / (GammaIdealFluid));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_dPdZ_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq) const {
  return ((1.0 - Vsq) * (GammaIdealFluid - 1.0) / GammaIdealFluid);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_dPdVsq_funcZVsq(CCTK_REAL Z, CCTK_REAL Vsq,
                                 cons_vars &cv) const {
  return ((-Z + cv.dens / (2.0 * sqrt(1.0 - Vsq))) * (GammaIdealFluid - 1.0) /
          GammaIdealFluid);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_2DNoble_f0(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL Ssq,
                            CCTK_REAL Bsq, CCTK_REAL BiSi) const {
  return (Vsq * (Bsq + Z) * (Bsq + Z) -
          (BiSi * BiSi * (Bsq + 2.0 * Z)) / (Z * Z) - Ssq);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_2DNoble_f1(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL Bsq,
                            CCTK_REAL BiSi, cons_vars &cv) const {
  CCTK_REAL Press = get_Press_funcZVsq(Z, Vsq, cv);
  return cv.tau + cv.dens - (Bsq / 2.0) * (1 + Vsq) +
         (BiSi * BiSi) / (2.0 * (Z * Z)) - Z + Press;

  //	printf("dum, press, tau, dens, Bsq, Vsq, BiSi, Z: %f, %f, %f, %f, %f,
  //%f, %f, %f \n", dum, Press, cv.tau, cv.dens, Bsq, Vsq, BiSi, Z);
  //        return cv.tau + cv.dens - Bsq / 2.0 * (1 + Vsq) +
  //               BiSi * BiSi / 2.0 / (Z * Z) - Z + Press;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_2DNoble_df0dZ(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL Bsq,
                               CCTK_REAL BiSi) const {
  return (2.0 * Vsq * (Bsq + Z) - 2.0 * BiSi * BiSi / (Z * Z) +
          2.0 * BiSi * BiSi * (Bsq + 2.0 * Z) / (Z * Z * Z));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_2DNoble_df0dVsq(CCTK_REAL Z, CCTK_REAL Bsq) const {
  return (Bsq + Z) * (Bsq + Z);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_2DNoble_df1dZ(CCTK_REAL Z, CCTK_REAL Vsq,
                               CCTK_REAL BiSi) const {
  return (-BiSi * BiSi / (Z * Z * Z) - 1.0 + get_dPdZ_funcZVsq(Z, Vsq));
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_2DNoble::get_2DNoble_df1dVsq(CCTK_REAL Z, CCTK_REAL Vsq, CCTK_REAL Bsq,
                                 cons_vars &cv) const {
  return (-Bsq / 2.0 + get_dPdVsq_funcZVsq(Z, Vsq, cv));
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_2DNoble::WZ2Prim(CCTK_REAL Z_Sol, CCTK_REAL vsq_Sol, CCTK_REAL Bsq,
                     CCTK_REAL BiSi, EOSType &eos_th, prim_vars &pv,
                     cons_vars &cv, const smat<CCTK_REAL, 3> &gup) const {
  CCTK_REAL W_Sol = 1.0 / sqrt(1.0 - vsq_Sol);

  pv.rho = cv.dens / W_Sol;

  pv.vel(X) =
      (gup(X, X) * cv.mom(X) + gup(X, Y) * cv.mom(Y) + gup(X, Z) * cv.mom(Z)) /
      (Z_Sol + Bsq);
  pv.vel(X) += BiSi * cv.dBvec(X) / (Z_Sol * (Z_Sol + Bsq));

  pv.vel(Y) =
      (gup(X, Y) * cv.mom(X) + gup(Y, Y) * cv.mom(Y) + gup(Y, Z) * cv.mom(Z)) /
      (Z_Sol + Bsq);
  pv.vel(Y) += BiSi * cv.dBvec(Y) / (Z_Sol * (Z_Sol + Bsq));

  pv.vel(Z) =
      (gup(X, Z) * cv.mom(X) + gup(Y, Z) * cv.mom(Y) + gup(Z, Z) * cv.mom(Z)) /
      (Z_Sol + Bsq);
  pv.vel(Z) += BiSi * cv.dBvec(Z) / (Z_Sol * (Z_Sol + Bsq));

  pv.w_lor = W_Sol;

  pv.eps = (Z_Sol * (1. - vsq_Sol) / pv.rho - 1.0) / GammaIdealFluid;

  pv.Ye = cv.dYe / cv.dens;

  pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.Bvec = cv.dBvec;

  const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
  pv.E = calc_contraction(gup, Elow);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_2DNoble::solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds,
                   cons_vars &cv, const smat<CCTK_REAL, 3> &glo,
                   CCTK_INT &c2p_succeeded) const {

  /* Calculate inverse of 3-metric */
  const CCTK_REAL spatial_detg = calc_det(glo);
  const CCTK_REAL sqrt_detg = sqrt(spatial_detg);
  const smat<CCTK_REAL, 3> gup = calc_inv(glo, spatial_detg);

  /* Undensitize the conserved vars */
  cv.dens /= sqrt_detg;
  cv.tau /= sqrt_detg;
  cv.mom /= sqrt_detg;
  cv.dBvec /= sqrt_detg;
  cv.dYe /= sqrt_detg;

  // compute primitive B seed from conserved B of current time step for better
  // guess
  pv_seeds.Bvec = cv.dBvec;

  const CCTK_REAL Ssq = get_Ssq_Exact(cv.mom, gup);
  const CCTK_REAL Bsq = get_Bsq_Exact(pv_seeds.Bvec, glo);
  const CCTK_REAL BiSi = get_BiSi_Exact(pv_seeds.Bvec, cv.mom);
  const vec<CCTK_REAL, 2> w_bsq = get_WLorentz_bsq_Seeds(
      pv_seeds.Bvec, pv_seeds.vel, glo); // this also recomputes pv_seeds.w_lor
  pv_seeds.w_lor = w_bsq(0);
  // const CCTK_REAL bsq = w_bsq(1);

  /* update rho seed from cv and wlor */
  // rho consistent with cv.rho should be better guess than rho from last
  // timestep
  pv_seeds.rho = cv.dens / pv_seeds.w_lor;

  /* get pressure seed from updated pv_seeds.rho */
  pv_seeds.press = eos_th.press_from_valid_rho_eps_ye(
      pv_seeds.rho, pv_seeds.eps, pv_seeds.Ye);

  /* get Z seed */
  CCTK_REAL Z_Seed =
      get_Z_Seed(pv_seeds.rho, pv_seeds.eps, pv_seeds.press, pv_seeds.w_lor);

  /* initialize unknowns for c2p, Z and vsq: */
  CCTK_REAL x[2];
  CCTK_REAL x_old[2];
  x[0] = fabs(Z_Seed);
  x[1] = (-1.0 + pv_seeds.w_lor * pv_seeds.w_lor) /
         (pv_seeds.w_lor * pv_seeds.w_lor);
  //      printf("x[0], x[1]: %f, %f \n", x[0], x[1]);

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
  c2p_succeeded = false; // false
  CCTK_INT k;
  for (k = 1; k <= maxIterations; k++) {
    fvec[0] = get_2DNoble_f0(x[0], x[1], Ssq, Bsq, BiSi);
    fvec[1] = get_2DNoble_f1(x[0], x[1], Bsq, BiSi, cv);
    fjac[0][0] = get_2DNoble_df0dZ(x[0], x[1], Bsq, BiSi);
    fjac[0][1] = get_2DNoble_df0dVsq(x[0], x[1]);
    fjac[1][0] = get_2DNoble_df1dZ(x[0], x[1], BiSi);
    fjac[1][1] = get_2DNoble_df1dVsq(x[0], x[1], Bsq, cv);
    detjac_inv = 1.0 / (fjac[0][0] * fjac[1][1] - fjac[0][1] * fjac[1][0]);
    dx[0] = -detjac_inv * (fjac[1][1] * fvec[0] - fjac[0][1] * fvec[1]);
    dx[1] = -detjac_inv * (-fjac[1][0] * fvec[0] + fjac[0][0] * fvec[1]);

    errf = 0.0;
    for (CCTK_INT i = 0; i < n; i++) {
      errf += fabs(fvec[i]);
      // printf("i, fabs(fvec[i]): %i, %f \n", i, fabs(fvec[i]));
    }

    /* save old values before calculating the new */
    x_old[0] = x[0];
    x_old[1] = x[1];

    for (CCTK_INT i = 0; i < n; i++) {
      x[i] += dx[i];
    }

    //	errf  = (x[0]==0.) ?  fabs(dx[0]) : fabs(dx[0]/x[0]);
    if (errf <= tolerance) {
      c2p_succeeded = true; // false
      //     printf("Noble c2p failed, c2p_succeeded set to 0. \n");
      //     printf("errf, tolerance: %f, %f \n", errf, tolerance);
      break;
    }

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
  }

  /* Calculate primitives from Z and W */
  CCTK_REAL Z_Sol = x[0];
  CCTK_REAL vsq_Sol = x[1];

  /* Write prims if C2P succeeded */
  /* TODO: report error */
  WZ2Prim(Z_Sol, vsq_Sol, Bsq, BiSi, eos_th, pv, cv, gup);
}

/* Destructor */
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_2DNoble::~c2p_2DNoble() {
  // How to destruct properly a vector?
}
} // namespace Con2PrimFactory

#endif
