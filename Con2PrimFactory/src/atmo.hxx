/*! \file atmo.hxx
\brief Class definition representing artificial atmosphere.
*/

#ifndef ATMO_HXX
#define ATMO_HXX

#include "prims.hxx"
#include "cons.hxx"

namespace Con2PrimFactory {

/// Class representing an artificial atmosphere.
struct atmosphere {
  const CCTK_REAL rho_atmo;
  const CCTK_REAL eps_atmo;
  const CCTK_REAL ye_atmo;
  const CCTK_REAL press_atmo;
  const CCTK_REAL rho_cut;

  atmosphere(const atmosphere &) = default;

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere(
      CCTK_REAL rho_, CCTK_REAL eps_, CCTK_REAL Ye_, CCTK_REAL press_,
      CCTK_REAL rho_cut_)
      : rho_atmo(rho_), eps_atmo(eps_), ye_atmo(Ye_), press_atmo(press_),
        rho_cut(rho_cut_) {}

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set(prim_vars &pv) const {

    pv.rho = rho_atmo;
    pv.eps = eps_atmo;
    pv.Ye = ye_atmo;
    pv.press = press_atmo;
    pv.vel(0) = 0.0;
    pv.vel(1) = 0.0;
    pv.vel(2) = 0.0;
    pv.w_lor = 1.0;
    pv.E(0) = 0.0;
    pv.E(1) = 0.0;
    pv.E(2) = 0.0;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set(prim_vars &pv, cons_vars &cv, const smat<CCTK_REAL, 3> &g) const {

    set(pv);
    const CCTK_REAL sqrt_detg = sqrt(calc_det(g));
    cv.dens = sqrt_detg * rho_atmo;
    cv.mom(0) = 0.0;
    cv.mom(1) = 0.0;
    cv.mom(2) = 0.0;
    cv.dYe = cv.dens * ye_atmo;
    const vec<CCTK_REAL, 3> &B_up = pv.Bvec;
    const vec<CCTK_REAL, 3> B_low = calc_contraction(g, B_up);
    CCTK_REAL Bsq = calc_contraction(B_up, B_low);
    cv.tau = cv.dens * eps_atmo + 0.5 * sqrt_detg * Bsq;
  }
};

} // namespace Con2PrimFactory
#endif
