#ifndef THERMO_STATE_HXX
#define THERMO_STATE_HXX

#include <cctk.h>

#include <cmath>

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

} // namespace EOSX

#endif
