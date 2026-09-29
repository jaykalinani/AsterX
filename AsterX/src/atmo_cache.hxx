#ifndef ASTERX_ATMO_CACHE_HXX
#define ASTERX_ATMO_CACHE_HXX

#include "atmo.hxx"

namespace AsterX {

// Filled once during global EOS validation, before grid kernels run. Loads
// are read-only: no host access to the managed EOS table is needed in the
// evolution loop, and concurrent grid calls do not mutate shared state.
struct atmo_cache {
  Con2PrimFactory::atmosphere atmo{};
  const void *cold = nullptr;
  const void *eos = nullptr;
  CCTK_REAL rho = 0.0;
  CCTK_REAL press = 0.0;
  CCTK_REAL temp = 0.0;
  CCTK_REAL Ye = 0.0;
  bool thermal = false;
  bool use_press = false;
  bool valid = false;

  template <typename EOSIDType, typename EOSType>
  void store(const EOSIDType *eos_1p, const EOSType *eos_3p,
              const Con2PrimFactory::atmosphere &state,
              const CCTK_REAL rho_abs_min, const CCTK_REAL p_atmo,
              const CCTK_REAL t_atmo, const CCTK_REAL Ye_atmo,
              const bool thermal_eos_atmo, const bool use_press_atmo) {
    atmo = state;
    cold = eos_1p;
    eos = eos_3p;
    rho = rho_abs_min;
    press = p_atmo;
    temp = t_atmo;
    Ye = Ye_atmo;
    thermal = thermal_eos_atmo;
    use_press = use_press_atmo;
    valid = true;
  }

  template <typename EOSIDType, typename EOSType>
  bool load(const EOSIDType *eos_1p, const EOSType *eos_3p,
             const CCTK_REAL rho_abs_min, const CCTK_REAL p_atmo,
             const CCTK_REAL t_atmo, const CCTK_REAL Ye_atmo,
             const CCTK_REAL n_rho_atmo, const CCTK_REAL n_press_atmo,
             const CCTK_REAL n_temp_atmo, const CCTK_REAL atmo_tol,
             const bool thermal_eos_atmo, const bool use_press_atmo,
             Con2PrimFactory::atmosphere &state) const {
    // Only the active grading exponents matter. Cold matching depends on
    // rho alone; its temperature and pressure exponents remain inactive.
    const bool graded = n_rho_atmo != 0.0 ||
        (thermal_eos_atmo &&
         (use_press_atmo ? n_press_atmo : n_temp_atmo) != 0.0);
    if (!valid || graded || cold != eos_1p || eos != eos_3p ||
        rho != rho_abs_min || press != p_atmo || temp != t_atmo ||
        Ye != Ye_atmo || thermal != thermal_eos_atmo ||
        use_press != use_press_atmo)
      return false;

    state = atmo;
    // atmo_tol is steerable and affects the cutoff, not thermodynamics.
    state.rho_cut = state.rho_atmo * (1.0 + atmo_tol);
    return true;
  }
};

inline atmo_cache cached_atmo;

} // namespace AsterX

#endif
