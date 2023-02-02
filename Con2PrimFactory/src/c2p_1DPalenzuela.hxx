#ifndef C2P_1DPALENZUELA_HXX
#define C2P_1DPALENZUELA_HXX

#include "c2p.hxx"
namespace Con2PrimFactory {

class c2p_1DPalenzuela : public c2p {
public:
  /* Some attributes */
  CCTK_REAL GammaIdealFluid;

  /* Constructor */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DPalenzuela(
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

  /* Called by 1DPalenzuela */
  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  xPalenzuelaToPrim(CCTK_REAL xPalenzuela_Sol, CCTK_REAL Ssq, CCTK_REAL Bsq,
                    CCTK_REAL BiSi, EOSType &eos_th, prim_vars &pv,
                    cons_vars &cv, const smat<CCTK_REAL, 3> &gup,
                    const smat<CCTK_REAL, 3> &glo) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE
      CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::pair<CCTK_REAL, CCTK_REAL>
      brent(CCTK_REAL Ssq, CCTK_REAL Bsq, CCTK_REAL BiSi, CCTK_INT min_bits,
            CCTK_INT max_iters, CCTK_INT &iters, EOSType &eos_th,
            cons_vars &cv) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  funcRoot_1DPalenzuela(CCTK_REAL Ssq, CCTK_REAL Bsq, CCTK_REAL BiSi,
                        CCTK_REAL x, EOSType &eos_th, cons_vars &cv) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds, cons_vars &cv,
        const smat<CCTK_REAL, 3> &glo, CCTK_INT &c2p_succeeded) const;

  /* Destructor */
  CCTK_HOST CCTK_DEVICE ~c2p_1DPalenzuela();
};

/* Constructor */
template <typename EOSType>
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DPalenzuela::c2p_1DPalenzuela(
        EOSType &eos_th, CCTK_INT maxIter, CCTK_REAL tol) {

  GammaIdealFluid = eos_th.gamma;
  maxIterations = maxIter;
  tolerance = tol;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DPalenzuela::get_Ssq_Exact(vec<CCTK_REAL, 3> &mom,
                                const smat<CCTK_REAL, 3> &gup) const {

  CCTK_REAL Ssq;
  /* calculate S_squared */
  Ssq = mom(X) * (gup(X, X) * mom(X) + gup(X, Y) * mom(Y) + gup(X, Z) * mom(Z));
  Ssq +=
      mom(Y) * (gup(X, Y) * mom(X) + gup(Y, Y) * mom(Y) + gup(Y, Z) * mom(Z));
  Ssq +=
      mom(Z) * (gup(X, Z) * mom(X) + gup(Y, Z) * mom(Y) + gup(Z, Z) * mom(Z));

  if ((Ssq < 0.) && (fabs(Ssq) < 1.0e-13)) {
    Ssq = fabs(Ssq);
  }
  return Ssq;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DPalenzuela::get_Bsq_Exact(vec<CCTK_REAL, 3> &B_up,
                                const smat<CCTK_REAL, 3> &glo) const {

  vec<CCTK_REAL, 3> B_low = calc_contraction(glo, B_up);
  return calc_contraction(B_low, B_up); // Bsq
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DPalenzuela::get_BiSi_Exact(vec<CCTK_REAL, 3> &Bvec,
                                 vec<CCTK_REAL, 3> &mom) const {

  return Bvec(X) * mom(X) + Bvec(Y) * mom(Y) + Bvec(Z) * mom(Z); // BiSi
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 2>
c2p_1DPalenzuela::get_WLorentz_bsq_Seeds(vec<CCTK_REAL, 3> &B_up,
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

/* Called by 1DPalenzuela */

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DPalenzuela::xPalenzuelaToPrim(CCTK_REAL xPalenzuela_Sol, CCTK_REAL Ssq,
                                    CCTK_REAL Bsq, CCTK_REAL BiSi,
                                    EOSType &eos_th, prim_vars &pv,
                                    cons_vars &cv,
                                    const smat<CCTK_REAL, 3> &gup,
                                    const smat<CCTK_REAL, 3> &glo) const {
  const CCTK_REAL qPalenzuela = cv.tau / cv.dens;
  const CCTK_REAL rPalenzuela = Ssq / pow(cv.dens, 2);
  const CCTK_REAL sPalenzuela = Bsq / cv.dens;
  const CCTK_REAL tPalenzuela = BiSi / pow(cv.dens, 3. / 2.);

  // (i)
  CCTK_REAL Wminus2 =
      1.0 -
      (xPalenzuela_Sol * xPalenzuela_Sol * rPalenzuela +
       (2.0 * xPalenzuela_Sol + sPalenzuela) * tPalenzuela * tPalenzuela) /
          (xPalenzuela_Sol * xPalenzuela_Sol * (xPalenzuela_Sol + sPalenzuela) *
           (xPalenzuela_Sol + sPalenzuela));
  Wminus2 = fmin(fmax(Wminus2, 1e-10), 1 - 1e-10);
  const CCTK_REAL W_sol = pow(Wminus2, -0.5);

  // (ii)
  pv.rho = cv.dens / W_sol;

  // (iii)
  pv.eps = W_sol - 1.0 + (1.0 - W_sol * W_sol) * xPalenzuela_Sol / W_sol +
           W_sol * (qPalenzuela - sPalenzuela +
                    tPalenzuela * tPalenzuela /
                        (2 * xPalenzuela_Sol * xPalenzuela_Sol) +
                    sPalenzuela / (2.0 * W_sol * W_sol));

  // (iv)
  // CCTK_REAL P_loc = get_Press_funcRhoEps(rho_loc, eps_loc);

  // Taken from WZ2Prim (2DNRNoble)
  CCTK_REAL Z_Sol = xPalenzuela_Sol * pv.rho * W_sol;

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

  pv.w_lor = W_sol;

  vec<CCTK_REAL, 3> v_low = calc_contraction(glo, pv.vel);
  CCTK_REAL vsq_Sol = calc_contraction(v_low, pv.vel);

  pv.eps = (Z_Sol * (1. - vsq_Sol) / pv.rho - 1.0) / GammaIdealFluid;

  pv.Ye = cv.dYe / cv.dens;

  pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

  pv.Bvec = cv.dBvec;

  const vec<CCTK_REAL, 3> Elow = calc_cross_product(pv.Bvec, pv.vel);
  pv.E = calc_contraction(gup, Elow);
}

template <typename EOSType>
CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline std::pair<CCTK_REAL, CCTK_REAL>
c2p_1DPalenzuela::brent(CCTK_REAL Ssq, CCTK_REAL Bsq, CCTK_REAL BiSi,
                        CCTK_INT min_bits, CCTK_INT max_iters, CCTK_INT &iters,
                        EOSType &eos_th, cons_vars &cv) const {
  using std::abs, std::min, std::max, std::swap;

  CCTK_REAL qPalenzuela = cv.tau / cv.dens;
  CCTK_REAL sPalenzuela = Bsq / cv.dens;

  CCTK_REAL xPalenzuela_lowerBound = 1.0 + qPalenzuela - sPalenzuela;
  CCTK_REAL xPalenzuela_upperBound = 2.0 + 2.0 * qPalenzuela - sPalenzuela;

  CCTK_REAL a = xPalenzuela_lowerBound;
  CCTK_REAL b = xPalenzuela_upperBound;

  auto tol = boost::math::tools::eps_tolerance<CCTK_REAL>(min_bits);

  iters = 0;
  auto fa = funcRoot_1DPalenzuela(Ssq, Bsq, BiSi, a, eos_th, cv);
  auto fb = funcRoot_1DPalenzuela(Ssq, Bsq, BiSi, b, eos_th, cv);
  if (abs(fa) < abs(fb)) {
    swap(a, b);
    swap(fa, fb);
  }
  if (fb == 0)
    return {b, b};
  if (fa * fb >= 0) {
    // Root is not bracketed
    iters = max_iters;
    return {min(a, b), max(a, b)};
  }
  CCTK_REAL c = a;
  auto fc = fa;
  bool mflag = true;
  CCTK_REAL d{};
  while (fb != 0 && !tol(a, b) && iters < max_iters) {
    CCTK_REAL s;
    if (fa != fc && fb != fc)
      // inverse quadratic interpolation
      s = (a * fb * fc) / ((fa - fb) * (fa - fc)) +
          (b * fa * fc) / ((fb - fa) * (fb - fc)) +
          (c * fa * fb) / ((fc - fa) * (fc - fb));
    else
      // secant method
      s = (a + b) / 2 - (fa + fb) / 2 * (b - a) / (fb - fa);

    CCTK_REAL u = (3 * a + b) / 4;
    CCTK_REAL v = b;

    if (u > v)
      swap(u, v);

    bool cond1 = !(u <= s && s <= v);
    bool cond2 = mflag && abs(s - b) >= abs(b - c) / 2;
    bool cond3 = !mflag && abs(s - b) >= abs(c - d) / 2;
    bool cond4 = mflag && tol(c, b);
    bool cond5 = !mflag && tol(c, d);
    if (cond1 || cond2 || cond3 || cond4 || cond5) {
      // bisection
      s = (a + b) / 2;
      mflag = true;
    } else {
      mflag = false;
    }
    auto fs = funcRoot_1DPalenzuela(Ssq, Bsq, BiSi, s, eos_th, cv);
    // `d` is assigned for the first time here; it won't be used above on the
    // first iteration because `mflag` is set
    d = c;
    c = b;
    fc = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }
    // CCTK_VINFO("iters=%d mflag=%d   a=%.17g b=%.17g c=%.17g d=%.17g fa=%.17g"
    //            "fb=%.17g fc=%.17g",
    //            iters, int(mflag), double(a), double(b), double(c), double(d),
    //            double(fa), double(fb), double(fc));
    if (fa * fb >= 0) {
      return {min(a, b), max(a, b)};
    }
    if (abs(fa) < abs(fb)) {
      swap(a, b);
      swap(fa, fb);
    }
    ++iters;
  }

  if (fb == 0)
    return {b, b};
  return {min(a, b), max(a, b)};
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
c2p_1DPalenzuela::funcRoot_1DPalenzuela(CCTK_REAL Ssq, CCTK_REAL Bsq,
                                        CCTK_REAL BiSi, CCTK_REAL x,
                                        EOSType &eos_th, cons_vars &cv) const {
  // computes f(x) from x and q,r,s,t
  const CCTK_REAL qPalenzuela = cv.tau / cv.dens;
  const CCTK_REAL rPalenzuela = Ssq / pow(cv.dens, 2);
  const CCTK_REAL sPalenzuela = Bsq / cv.dens;
  const CCTK_REAL tPalenzuela = BiSi / pow(cv.dens, 3. / 2.);

  // (i)
  CCTK_REAL Wminus2 =
      1.0 - (x * x * rPalenzuela +
             (2.0 * x + sPalenzuela) * tPalenzuela * tPalenzuela) /
                (x * x * (x + sPalenzuela) * (x + sPalenzuela));
  Wminus2 = fmin(fmax(Wminus2, 1e-10), 1 - 1e-10);
  const CCTK_REAL W_loc = pow(Wminus2, -0.5);

  // (ii)
  CCTK_REAL rho_loc = cv.dens / W_loc;
  CCTK_REAL Ye_loc = cv.dYe / cv.dens;

  // (iii)
  CCTK_REAL eps_loc = W_loc - 1.0 + (1.0 - W_loc * W_loc) * x / W_loc +
                      W_loc * (qPalenzuela - sPalenzuela +
                               tPalenzuela * tPalenzuela / (2 * x * x) +
                               sPalenzuela / (2 * W_loc * W_loc));

  // (iv)
  CCTK_REAL P_loc =
      eos_th.press_from_valid_rho_eps_ye(rho_loc, eps_loc, Ye_loc);

  return (x - (1.0 + eps_loc + P_loc / rho_loc) * W_loc);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p_1DPalenzuela::solve(EOSType &eos_th, prim_vars &pv, prim_vars &pv_seeds,
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

  // Send con2primFactory object as reference to modify it,
  // and because we can not instantiate abstract class

  c2p_succeeded = false; // 0 //false

  // Find x, this is the recovery process
  const CCTK_INT minbits = std::numeric_limits<CCTK_REAL>::digits - 4;
  const CCTK_INT maxiters = maxIterations;
  CCTK_INT iters;
  // CCTK_REAL froot = funcRoot_1DPalenzuela();
  std::pair<CCTK_REAL, CCTK_REAL> result =
      brent(Ssq, Bsq, BiSi, minbits, maxiters, iters, eos_th, cv);

  // Pick best solution
  CCTK_REAL xPalenzuela_Sol;
  if (abs(funcRoot_1DPalenzuela(Ssq, Bsq, BiSi, result.first, eos_th, cv)) <
      abs(funcRoot_1DPalenzuela(Ssq, Bsq, BiSi, result.second, eos_th, cv))) {
    xPalenzuela_Sol = result.first;
  } else {
    xPalenzuela_Sol = result.second;
  }

  // Check solution and calculate primitives
  // TODO:check if to pass result.first or xPalenzuela_Sol
  if (iters < maxiters &&
      abs(funcRoot_1DPalenzuela(Ssq, Bsq, BiSi, result.first, eos_th, cv)) <
          tolerance) {
    //        printf("Palenzuela C2P failed, c2p_succeeded set to 0. \n");
    c2p_succeeded = true; // true
  }

  xPalenzuelaToPrim(xPalenzuela_Sol, Ssq, Bsq, BiSi, eos_th, pv, cv, gup, glo);
  return;
}

/* Destructor */
CCTK_HOST CCTK_DEVICE
    CCTK_ATTRIBUTE_ALWAYS_INLINE inline c2p_1DPalenzuela::~c2p_1DPalenzuela() {
  // How to destruct properly a vector?
}
} // namespace Con2PrimFactory

#endif
