#ifndef C2P_1DENTROPY_HXX
#define C2P_1DENTROPY_HXX

#include "c2p.hxx"
#include "roots.hxx"
namespace Con2PrimFactory {

class c2p_1DEntropy : public c2p {
public:
  /* Some attributes */
  CCTK_REAL GammaIdealFluid;

  /* Constructor */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DEntropy(
      EOSType &eos_th, atmosphere &atm, CCTK_INT maxIter, CCTK_REAL tol,
      CCTK_REAL rho_str, CCTK_REAL vwlim, CCTK_REAL B_lim, bool ye_len);

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Ssq_Exact(vec<CCTK_REAL, 3> &mom, const smat<CCTK_REAL, 3> &gup) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up, const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec, vec<CCTK_REAL, 3> &mom) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
  get_WLorentz_vsq_bsq_Seeds(vec<CCTK_REAL, 3> &B_up, vec<CCTK_REAL, 3> &v_up,
                             const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set_to_nan(prim_vars &pv, cons_vars &cv) const;
  /* Called by 1DEntropy */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  xEntropyToPrim(CCTK_REAL xEntropy_Sol, CCTK_REAL Ssq, CCTK_REAL Bsq,
                 CCTK_REAL BiSi, EOSType &eos_th, prim_vars &pv, cons_vars cv,
                 const smat<CCTK_REAL, 3> &gup,
                 const smat<CCTK_REAL, 3> &glo) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  funcRoot_1DEntropy(CCTK_REAL Ssq, CCTK_REAL Bsq, CCTK_REAL BiSi, CCTK_REAL x,
                     EOSType &eos_th, cons_vars &cv) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds, cons_vars cv,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const;

  /* Destructor */
  CCTK_HOST CCTK_DEVICE ~c2p_1DEntropy();
};

/* Constructor */
template <typename EOSType>
CCTK_HOST CCTK_DEVICE
CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DEntropy::c2p_1DEntropy(
    EOSType &eos_th, atmosphere &atm, CCTK_INT maxIter, CCTK_REAL tol,
    CCTK_REAL rho_str, CCTK_REAL vwlim, CCTK_REAL B_lim, bool ye_len) {

  GammaIdealFluid = eos_th.gamma;
  maxIterations = maxIter;
  tolerance = tol;
  rho_strict = rho_str;
  ye_lenient = ye_len;
  vw_lim = vwlim;
  w_lim = sqrt(1.0 + vw_lim * vw_lim);
  v_lim = vw_lim / w_lim;
  Bsq_lim = B_lim * B_lim;
  atmo = atm;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DEntropy::get_Ssq_Exact(vec<CCTK_REAL, 3> &mom,
                             const smat<CCTK_REAL, 3> &gup) const {

  CCTK_REAL Ssq;
  /* calculate S_squared */
  Ssq = mom(X) * (gup(X, X) * mom(X) + gup(X, Y) * mom(Y) + gup(X, Z) * mom(Z));
  Ssq +=
      mom(Y) * (gup(X, Y) * mom(X) + gup(Y, Y) * mom(Y) + gup(Y, Z) * mom(Z));
  Ssq +=
      mom(Z) * (gup(X, Z) * mom(X) + gup(Y, Z) * mom(Y) + gup(Z, Z) * mom(Z));

  // if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13)) {
  //  Ssq = fabs(Ssq);
  //}
  return Ssq;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DEntropy::get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up,
                             const smat<CCTK_REAL, 3> &glo) const {

  vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);
  return calc_contraction(B_low, B_up); // Bsq
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DEntropy::get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec,
                              vec<CCTK_REAL, 3> &mom) const {

  return Bvec(X) * mom(X) + Bvec(Y) * mom(Y) + Bvec(Z) * mom(Z); // BiSi
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 3>
c2p_1DEntropy::get_WLorentz_vsq_bsq_Seeds(vec<CCTK_REAL, 3> &B_up,
                                          vec<CCTK_REAL, 3> &v_up,
                                          const smat<CCTK_REAL, 3> &glo) const {
  vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);
  CCTK_REAL vsq = calc_contraction(v_low, v_up);
  CCTK_REAL VdotB = calc_contraction(v_low, B_up);
  CCTK_REAL VdotBsq = VdotB * VdotB;
  CCTK_REAL Bsq = get_Bsq_Exact(B_up, glo);

  // if ((vsq < 0.) && (fabs(vsq) < 1.0e-13)) {
  //  vsq = fabs(vsq);
  //}

  CCTK_REAL w_lor = 1. / sqrt(1. - vsq);
  CCTK_REAL bsq = ((Bsq) / (w_lor * w_lor)) + VdotBsq;
  vec<CCTK_REAL, 3> w_vsq_bsq{ w_lor, vsq, bsq };

  return w_vsq_bsq; //{w_lor, vsq, bsq}
}

/* Called by 1DEntropy */

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DEntropy::xEntropyToPrim(CCTK_REAL xEntropy_Sol, CCTK_REAL Ssq,
                              CCTK_REAL Bsq, CCTK_REAL BiSi, EOSType &eos_th,
                              prim_vars &pv, cons_vars cv,
                              const smat<CCTK_REAL, 3> &gup,
                              const smat<CCTK_REAL, 3> &glo) const {
  // Density, entropy, Ye
  pv.rho = xEntropy_Sol;
  pv.entropy = cv.dS / cv.dens;
  pv.Ye = cv.dYe / cv.dens;

  // Lorentz factor
  pv.w_lor = cv.dens / xEntropy_Sol;

  // Pressure and epsilon
  pv.press =
      eos_th.press_from_valid_rho_entropy_ye(xEntropySol, pv.entropy, pv.Ye);
  pv.eps = eos_th.eps_from_valid_rho_entropy_ye(xEntropySol, pv.entropy, pv.Ye);

  // Taken from WZ2Prim (2DNRNoble)
  // Z_Sol = rho * h * w_lor * w_lor
  CCTK_REAL Z_Sol =
      (xEntropy_Sol + pv.eps * xEntropy_Sol + pv.press) * pv.w_lor * pv.w_lor;

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

  pv.Bvec = cv.dBvec;

  const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
  pv.E = calc_contraction(gup, Elow);
}

// See Appendix A.4 of https://arxiv.org/pdf/1112.0568
template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DEntropy::funcRoot_1DEntropy(CCTK_REAL Ssq, CCTK_REAL Bsq, CCTK_REAL BiSi,
                                  CCTK_REAL x, EOSType &eos_th,
                                  cons_vars &cv) const {

  // We already divided dens, dS and
  // dYe by sqrt(gamma)

  const CCTK_REAL ent_loc = cv.dS / dens;
  const CCTK_REAL ye_loc = cv.dYe / dens;

  // Compute h using entropy
  const CCTK_REAL press_loc =
      eos_th.press_from_valid_rho_entropy_ye(x, ent_loc, ye_loc);

  const CCTK_REAL eps_loc =
      eos_th.eps_from_valid_rho_entropy_ye(x, ent_loc, ye_loc);

  // Compute (A60) using
  // W = rho*h*lorentz*lorentz
  const CCTK_REAL lor_loc = cv.dens / x;
  const CCTK_REAL h_loc = (1.0 + eps_loc + press_loc / x);
  const CCTK_REAL W = x * h_loc * lor_loc * lor_loc;

  // Compute (A61)
  const CCTK_REAL Sfsq =
      (W * W * Ssq + BiSi * BiSi * (Bsq + 2.0 * W)) / ((W + Bsq) * (W + Bsq));

  // Compute (A62)
  const CCTK_REAL rho =
      cv.dens / sqrt(1.0 + Sfsq / (cv.dens * cv.dens * h_loc * h_loc));

  return x - rho;
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DEntropy::solve(EOSType &eos_th, prim_vars &pv, cons_vars cv,
                     const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const {

  ROOTSTAT status = ROOTSTAT::SUCCESS;
  rep.iters = 0;
  rep.adjust_cons = false;
  rep.set_atmo = false;
  rep.status = c2p_report::SUCCESS;

  /* Check validity of the 3-metric and compute its inverse */
  const CCTK_REAL spatial_detg = calc_det(glo);
  const CCTK_REAL sqrt_detg = sqrt(spatial_detg);
  // if ((!isfinite(sqrt_detg)) || (sqrt_detg <= 0)) {
  //  rep.set_invalid_detg(sqrt_detg);
  //  set_to_nan(pv, cv);
  //  return;
  //}

  // Check positive definiteness of spatial metric
  // Sylvester's criterion, see
  // https://en.wikipedia.org/wiki/Sylvester%27s_criterion
  const bool minor1{ glo(X, X) > 0.0 };
  const bool minor2{ glo(X, X) * glo(Y, Y) - glo(X, Y) * glo(X, Y) > 0.0 };
  const bool minor3{ spatial_detg > 0.0 };

  if (minor1 && minor2 && minor3) {
    continue;
  } else {
    rep.set_invalid_detg(sqrt_detg);
    set_to_nan(pv, cv);
    return;
  }

  const smat<CCTK_REAL, 3> gup = calc_inv(glo, spatial_detg);

  /* Undensitize the conserved vars */
  cv.dens /= sqrt_detg;
  // cv.tau /= sqrt_detg;
  cv.mom /= sqrt_detg;
  cv.dBvec /= sqrt_detg;
  cv.dYe /= sqrt_detg;
  cv.dS /= sqrt_detg;

  if (cv.dens <= atmo.rho_cut) {
    rep.set_atmo_set();
    atmo.set(pv, cv, glo);
    return;
  }

  // compute primitive B seed from conserved B of current time step for better
  // guess
  // pv_seeds.Bvec = cv.dBvec;

  const CCTK_REAL Ssq = get_Ssq_Exact(cv.mom, gup);
  const CCTK_REAL Bsq = get_Bsq_Exact(cv.dBvec, glo);
  const CCTK_REAL BiSi = get_BiSi_Exact(cv.dBvec, cv.mom);
  // const vec<CCTK_REAL, 3> w_vsq_bsq = get_WLorentz_vsq_bsq_Seeds(
  //    pv_seeds.Bvec, pv_seeds.vel, glo); // this also recomputes
  // pv_seeds.w_lor
  // pv_seeds.w_lor = w_vsq_bsq(0);
  // CCTK_REAL vsq_seed = w_vsq_bsq(1);
  // const CCTK_REAL bsq = w_vsq_bsq(2);

  // CHECK
  // if ((!isfinite(cv.dens)) || (!isfinite(Ssq)) || (!isfinite(Bsq)) ||
  //    (!isfinite(BiSi)) || (!isfinite(cv.dYe))) {
  //  rep.set_nans_in_cons(cv.dens, Ssq, Bsq, BiSi, cv.dYe);
  //  set_to_nan(pv, cv);
  //  return;
  //}

  // if (Bsq < 0) {
  //  rep.set_neg_Bsq(Bsq);
  //  set_to_nan(pv, cv);
  //  return;
  //}

  // CHECK
  if (Bsq > Bsq_lim) {
    rep.set_B_limit(Bsq);
    set_to_nan(pv, cv);
    return;
  }

  /* update rho seed from cv and wlor */
  // rho consistent with cv.rho should be better guess than rho from last
  // timestep
  // pv_seeds.rho = cv.dens / pv_seeds.w_lor;

  // See Appendix A.4 of https://arxiv.org/pdf/1112.0568
  // Compute (A59)
  CCTK_REAL a = cv.dens / sqrt(1.0 + Ssq / (cv.dens * cv.dens));
  CCTK_REAL b = cv.dens;
  auto fn = [&](auto x) {
    return funcRoot_1DEntropy(Ssq, Bsq, BiSi, x, eos_th, cv);
  };

  // Important!
  // Algo::brent terminates if the following accuracy is achieved
  // abs(x - y) <= eps * min(abs(x), abs(y)),
  // where x and y are the values of the bracket and
  // eps = std::ldexp(1, -minbits) = 1 * 2^{-minbits}
  // This should probably be changed in Algo::brent

  // We want to set the tolerance to its correct parameter
  constexpr CCTK_REAL log2 = std::log(2.0);
  const CCTK_INT minbits = int(abs(std::log(tolerance)) / log2);
  const CCTK_REAL tolerance_0 = std::ldexp(double(1.0), -minbits);
  // Old code:
  // const CCTK_INT minbits = std::numeric_limits<CCTK_REAL>::digits - 4;
  const CCTK_INT maxiters = maxIterations;

  auto result = Algo::brent(fn, a, b, minbits, maxiters, rep.iters);

  CCTK_REAL xEntropy_Sol = 0.5 * (result.first + result.second);

  // Check solution and calculate primitives
  //  if (rep.iters < maxiters && abs(fn(xEntropy_Sol)) < tolerance) {
  if (abs(result.first - result.second) <=
      tolerance_0 * min(abs(result.first), abs(result.second))) {
    rep.status = c2p_report::SUCCESS;
    status = ROOTSTAT::SUCCESS;
  } else {
    // set status to root not converged
    rep.set_root_conv();
    status = ROOTSTAT::NOT_CONVERGED;
  }

  xEntropyToPrim(xEntropy_Sol, Ssq, Bsq, BiSi, eos_th, pv, cv, gup, glo);

  // Lower velocity
  const vec<CCTK_REAL, 3> v_low = calc_contraction(glo, pv.vel);
  // Computing b^t : this is b^0 * alp
  const CCTK_REAL bst = pv.w_lor * calc_contraction(pv.Bvec, v_low);
  // Computing b^mu b_mu
  const CCTK_REAL bs2 = (Bsq + bst * bst) / (pv.w_lor * pv.w_lor);
  // Recompute tau
  cv.tau = (pv.w_lor * pv.w_lor * (pv.rho * (1.0 + pv.eps) + pv.press + bs2) -
            (pv.press + 0.5 * bs2) - bst * bst) -
           cv.dens;

  // ----------
  // Floor and ceiling for rho and velocity
  // Keeps pressure the same and changes eps
  // ----------

  // set to atmo if computed rho is below floor density
  if (pv.rho < atmo.rho_cut) {
    rep.set_atmo_set();
    atmo.set(pv, cv, glo);
    return;
  }

  // check if computed velocities are within the specified limit
  CCTK_REAL vsq_Sol = calc_contraction(v_low, pv.vel);
  CCTK_REAL sol_v = sqrt(vsq_Sol);
  if (sol_v > v_lim) {
    /*
    printf("(sol_v > v_lim) is true! \n");
    printf("sol_v, v_lim: %26.16e, %26.16e \n", sol_v, v_lim);
    */
    // add mass, keeps conserved density D
    pv.rho = cv.dens / w_lim;
    pv.eps = eos_th.eps_from_valid_rho_press_ye(pv.rho, pv.press, pv.Ye);
    pv.entropy =
        eos_th.entropy_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    // if (pv.rho >= rho_strict) {
    //  rep.set_speed_limit({ sol_v, sol_v, sol_v });
    //  set_to_nan(pv, cv);
    //  return;
    //}
    pv.vel *= v_lim / sol_v;
    pv.w_lor = w_lim;
    // pv.eps = std::min(std::max(eos_th.rgeps.min, pv.eps),
    // eos_th.rgeps.max);
    // pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  if (pv.rho > rho_strict) {
    rep.adjust_cons = true;
    // remove mass, changes conserved density D
    pv.rho = rho_strict;
    pv.eps = eos_th.eps_from_valid_rho_press_ye(rho_strict, pv.press, pv.Ye);
    pv.entropy =
        eos_th.entropy_from_valid_rho_eps_ye(rho_strict, pv.eps, pv.Ye);
  }

  // ----------
  // Floor and ceiling for eps
  // Keeps rho the same and changes press
  // ----------

  // check the validity of the computed eps
  auto rgeps = eos_th.range_eps_from_valid_rho_ye(pv.rho, pv.Ye);
  if (pv.eps > rgeps.max) {
    // printf("(pv.eps > rgeps.max) is true, adjusting cons.. \n");
    rep.adjust_cons = true;
    // if (pv.rho >= rho_strict) {
    //  rep.set_range_eps(pv.eps); // sets adjust_cons to false by default
    //  rep.adjust_cons = true;
    //  set_to_nan(pv, cv);
    //  return;
    //}
    pv.eps = rgeps.max;
    pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.entropy = eos_th.entropy_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
  } else if (pv.eps < rgeps.min) {
    /*
    printf(
        "(pv.eps < rgeps.min) is true! pv.eps, rgeps.min: %26.16e, %26.16e
    \n",
        pv.eps, rgeps.min);
    printf(" Not adjusting cons.. \n");
    */
    // rep.set_range_eps(rgeps.min); // sets adjust_cons to true
    rep.adjust_cons = true;
    pv.eps = rgeps.min;
    pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.entropy = eos_th.entropy_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
  }

  // TODO: check validity for Ye

  // Recompute cons if prims have been adjusted
  if (rep.adjust_cons) {
    cv.from_prim(pv, glo);
  }
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DPalenzuela::set_to_nan(prim_vars &pv, cons_vars &cv) const {
  pv.set_to_nan();
  cv.set_to_nan();
}

/* Destructor */
CCTK_HOST CCTK_DEVICE
CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DPalenzuela::~c2p_1DPalenzuela() {
  // How to destruct properly a vector?
}
} // namespace Con2PrimFactory

#endif
