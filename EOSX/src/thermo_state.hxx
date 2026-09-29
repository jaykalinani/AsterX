#ifndef THERMO_STATE_HXX
#define THERMO_STATE_HXX

#include <cctk.h>

#include <cmath>
#include <limits>

#include "eos_3p.hxx"

namespace EOSX {

struct thermo_state {
  CCTK_REAL rho;
  CCTK_REAL eps;
  CCTK_REAL Ye;
  CCTK_REAL press;
  CCTK_REAL temperature;
  CCTK_REAL kappa;
  CCTK_REAL cs2;
};

struct thermo_state_derivs {
  thermo_state state;
  CCTK_REAL dpdrho;
  CCTK_REAL dpdeps;
  bool enthalpy_clipped;
};

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
limit_to_range(const CCTK_REAL value, const eos_3p::range &range) {
  return fmin(fmax(value, range.min), range.max);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline thermo_state
state_from_rho_temp_ye(const EOSType *eos, const CCTK_REAL rho,
                       const CCTK_REAL temperature, const CCTK_REAL Ye) {
  thermo_state state;
  state.rho = limit_to_range(rho, eos->rgrho);
  state.temperature = limit_to_range(temperature, eos->rgtemp);
  state.Ye = limit_to_range(Ye, eos->rgye);
  state.eps =
      eos->eps_from_rho_temp_ye(state.rho, state.temperature, state.Ye);
  state.press =
      eos->press_from_rho_temp_ye(state.rho, state.temperature, state.Ye);
  state.kappa = eos->kappa_from_rho_eps_ye(state.rho, state.eps, state.Ye);
  const CCTK_REAL csound =
      eos->csnd_from_rho_temp_ye(state.rho, state.temperature, state.Ye);
  state.cs2 = csound * csound;
  return state;
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline thermo_state
state_from_rho_eps_ye(const EOSType *eos, const CCTK_REAL rho,
                      const CCTK_REAL eps, const CCTK_REAL Ye) {
  thermo_state state;
  state.rho = limit_to_range(rho, eos->rgrho);
  state.Ye = limit_to_range(Ye, eos->rgye);
  const auto eps_range = eos->range_eps_from_rho_ye(state.rho, state.Ye);
  state.eps = limit_to_range(eps, eps_range);
  state.press = eos->press_from_rho_eps_ye(state.rho, state.eps, state.Ye);
  state.temperature =
      eos->temp_from_rho_eps_ye(state.rho, state.eps, state.Ye);
  state.kappa = eos->kappa_from_rho_eps_ye(state.rho, state.eps, state.Ye);
  const CCTK_REAL csound =
      eos->csnd_from_rho_eps_ye(state.rho, state.eps, state.Ye);
  state.cs2 = csound * csound;
  return state;
}

// Only use pressure-primary closure for EOS implementations with a unique
// pressure inversion, such as the ideal-gas EOS.
template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline thermo_state
state_from_rho_press_ye(const EOSType *eos, const CCTK_REAL rho,
                        const CCTK_REAL press, const CCTK_REAL Ye) {
  const CCTK_REAL rho_limited = limit_to_range(rho, eos->rgrho);
  const CCTK_REAL Ye_limited = limit_to_range(Ye, eos->rgye);
  const CCTK_REAL eps =
      eos->eps_from_rho_press_ye(rho_limited, press, Ye_limited);
  return state_from_rho_eps_ye(eos, rho_limited, eps, Ye_limited);
}

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline thermo_state_derivs
state_from_rho_enthalpy_ye(const EOSType *eos, const CCTK_REAL rho,
                           const CCTK_REAL enthalpy, const CCTK_REAL Ye) {
  const CCTK_REAL rho_limited = limit_to_range(rho, eos->rgrho);
  const CCTK_REAL Ye_limited = limit_to_range(Ye, eos->rgye);
  const auto eps_range =
      eos->range_eps_from_rho_ye(rho_limited, Ye_limited);

  auto enthalpy_from_eps = [&](const CCTK_REAL eps) {
    CCTK_REAL epsL = eps;
    const CCTK_REAL press =
        eos->press_from_rho_eps_ye(rho_limited, epsL, Ye_limited);
    return 1.0 + epsL + press / rho_limited;
  };

  const CCTK_REAL hmin = enthalpy_from_eps(eps_range.min);
  const CCTK_REAL hmax = enthalpy_from_eps(eps_range.max);
  CCTK_REAL eps;
  bool enthalpy_clipped = false;
  if (enthalpy <= hmin) {
    eps = eps_range.min;
    enthalpy_clipped = enthalpy < hmin;
  } else if (enthalpy >= hmax) {
    eps = eps_range.max;
    enthalpy_clipped = enthalpy > hmax;
  } else {
    CCTK_REAL eps_lo = eps_range.min;
    CCTK_REAL eps_hi = eps_range.max;
    eps = eps_lo + (eps_hi - eps_lo) * (enthalpy - hmin) / (hmax - hmin);
    // Stable EOSs have dh/deps = 1 + (dP/deps)/rho > 0. Keep a bracket
    // around the solution and fall back to its midpoint if a Newton step
    // would leave the local EOS range.
    const CCTK_REAL htol =
        32.0 * std::numeric_limits<CCTK_REAL>::epsilon() *
        fmax(1.0, fabs(enthalpy));
    for (CCTK_INT n = 0; n < 32; ++n) {
      CCTK_REAL press;
      CCTK_REAL dpdrho;
      CCTK_REAL dpdeps;
      eos->press_derivs_from_rho_eps_ye(press, dpdrho, dpdeps, rho_limited,
                                        eps, Ye_limited);
      const CCTK_REAL f = 1.0 + eps + press / rho_limited - enthalpy;
      if (fabs(f) <= htol) {
        break;
      }
      if (f < 0.0) {
        eps_lo = eps;
      } else {
        eps_hi = eps;
      }

      const CCTK_REAL dhdeps = 1.0 + dpdeps / rho_limited;
      const CCTK_REAL eps_newton = eps - f / dhdeps;
      const bool use_newton = std::isfinite(eps_newton) &&
                              eps_newton > eps_lo && eps_newton < eps_hi;
      eps = use_newton ? eps_newton : 0.5 * (eps_lo + eps_hi);
    }
  }

  thermo_state_derivs result;
  result.state =
      state_from_rho_eps_ye(eos, rho_limited, eps, Ye_limited);
  eos->press_derivs_from_rho_eps_ye(
      result.state.press, result.dpdrho, result.dpdeps, result.state.rho,
      result.state.eps, result.state.Ye);
  result.enthalpy_clipped = enthalpy_clipped;
  return result;
}

} // namespace EOSX

#endif
