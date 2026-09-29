/*! \file atmo.hxx
\brief Class definition representing artificial atmosphere.
*/

#ifndef ATMO_HXX
#define ATMO_HXX

#include "prims.hxx"
#include "cons.hxx"

#include "aster_utils.hxx"
#include "thermo_state.hxx"

namespace Con2PrimFactory {

/// Class representing an artificial atmosphere.
struct atmosphere {
  CCTK_REAL rho_atmo;
  CCTK_REAL eps_atmo;
  CCTK_REAL ye_atmo;
  CCTK_REAL press_atmo;
  CCTK_REAL temp_atmo;
  CCTK_REAL entropy_atmo;
  CCTK_REAL rho_cut;

  CCTK_DEVICE
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere() = default;

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere(
      CCTK_REAL rho_, CCTK_REAL eps_, CCTK_REAL Ye_, CCTK_REAL press_,
      CCTK_REAL temp_, CCTK_REAL entropy_, CCTK_REAL rho_cut_)
      : rho_atmo(rho_), eps_atmo(eps_), ye_atmo(Ye_), press_atmo(press_),
        temp_atmo(temp_), entropy_atmo(entropy_), rho_cut(rho_cut_) {}

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere &
  operator=(const atmosphere &other) {
    if (this == &other)
      return *this; // Handle self-assignment
    // Copy data members from 'other' to 'this'
    rho_atmo = other.rho_atmo;
    eps_atmo = other.eps_atmo;
    ye_atmo = other.ye_atmo;
    press_atmo = other.press_atmo;
    temp_atmo = other.temp_atmo;
    entropy_atmo = other.entropy_atmo;
    rho_cut = other.rho_cut;
    return *this;
  }

  CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  set(prim_vars &pv) const {

    pv.rho = rho_atmo;
    pv.eps = eps_atmo;
    pv.Ye = ye_atmo;
    pv.press = press_atmo;
    pv.temperature = temp_atmo;
    pv.entropy = entropy_atmo;
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
    cv.DYe = cv.dens * ye_atmo;
    cv.DEnt = cv.dens * entropy_atmo;
    const vec<CCTK_REAL, 3> &B_up = pv.Bvec;
    const vec<CCTK_REAL, 3> B_low = calc_contraction(g, B_up);
    CCTK_REAL Bsq = calc_contraction(B_up, B_low);
    cv.tau = cv.dens * eps_atmo + 0.5 * sqrt_detg * Bsq;
  }
};

CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
graded_atmosphere_value(const CCTK_REAL value_at_inner_radius,
                        const CCTK_REAL radial_distance,
                        const CCTK_REAL inner_radius,
                        const CCTK_REAL exponent) {
  return radial_distance > inner_radius
             ? value_at_inner_radius *
                   pow(inner_radius / radial_distance, exponent)
             : value_at_inner_radius;
}

template <typename EOSIDType, typename EOSType>
CCTK_DEVICE CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline atmosphere
make_atmosphere(const EOSIDType *eos_1p, const EOSType *eos_3p,
                const CCTK_REAL radial_distance,
                const CCTK_REAL rho_abs_min, const CCTK_REAL p_atmo,
                const CCTK_REAL t_atmo, const CCTK_REAL Ye_atmo,
                const CCTK_REAL r_atmo, const CCTK_REAL n_rho_atmo,
                const CCTK_REAL n_press_atmo,
                const CCTK_REAL n_temp_atmo, const CCTK_REAL atmo_tol,
                const bool thermal_eos_atmo, const bool use_press_atmo) {
  const CCTK_REAL rho_requested = graded_atmosphere_value(
      rho_abs_min, radial_distance, r_atmo, n_rho_atmo);
  const CCTK_REAL rho_atm =
      EOSX::limit_to_range(rho_requested, eos_3p->rgrho);
  const CCTK_REAL Ye_atm = EOSX::limit_to_range(Ye_atmo, eos_3p->rgye);

  EOSX::thermo_state state;
  if (!thermal_eos_atmo) {
    // Match the atmosphere to the cold initial-data EOS, then close the
    // complete state with the evolution EOS. This preserves the traditional
    // polytropic ideal-gas atmosphere and its density-driven grading.
    const CCTK_REAL gm1 = eos_1p->gm1_from_rho(rho_atm);
    const CCTK_REAL eps_cold = eos_1p->sed_from_gm1(gm1);
    state = EOSX::state_from_rho_eps_ye(eos_3p, rho_atm, eps_cold, Ye_atm);
  } else if (use_press_atmo) {
    // Pressure-primary closure is valid only for EOS implementations with a
    // unique pressure inversion. Startup validation rejects unsupported uses.
    const CCTK_REAL press_requested = graded_atmosphere_value(
        p_atmo, radial_distance, r_atmo, n_press_atmo);
    const auto state_at_temp_min = EOSX::state_from_rho_temp_ye(
        eos_3p, rho_atm, eos_3p->rgtemp.min, Ye_atm);
    const CCTK_REAL press_atm =
        fmax(press_requested, state_at_temp_min.press);
    state = EOSX::state_from_rho_press_ye(eos_3p, rho_atm, press_atm, Ye_atm);
  } else {
    const CCTK_REAL temp_requested = graded_atmosphere_value(
        t_atmo, radial_distance, r_atmo, n_temp_atmo);
    state = EOSX::state_from_rho_temp_ye(eos_3p, rho_atm, temp_requested,
                                         Ye_atm);
  }

  const CCTK_REAL rho_cut = state.rho * (1.0 + atmo_tol);
  return atmosphere(state.rho, state.eps, state.Ye, state.press,
                    state.temperature, state.kappa, rho_cut);
}

} // namespace Con2PrimFactory
#endif
