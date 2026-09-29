#include <cctk.h>
#include <cctk_Arguments.h>

#include "setup_eos.hxx"
#include "thermo_state.hxx"

namespace AsterSeeds {
using namespace EOSX;

namespace {
template <typename EOSType> void test_closure(const EOSType *eos) {
  for (const CCTK_REAL f : {0.25, 0.5, 0.75}) {
    const CCTK_REAL rho = exp((1.0 - f) * log(eos->rgrho.min) +
                             f * log(eos->rgrho.max));
    const CCTK_REAL temp = eos->rgtemp.min + f *
        (eos->rgtemp.max - eos->rgtemp.min);
    const CCTK_REAL Ye = eos->rgye.min + f * (eos->rgye.max - eos->rgye.min);
    // The beta-floor adjustment changes rho at fixed T and Ye. Close from
    // that same authority, including kappa rather than physical entropy.
    const auto state = state_from_rho_temp_ye(eos, rho, temp, Ye);
    const auto changed = state_from_rho_temp_ye(eos, 1.1 * rho, temp, Ye);
    CCTK_REAL eps = changed.eps;
    const CCTK_REAL kappa = eos->kappa_from_rho_eps_ye(
        changed.rho, eps, changed.Ye);
    const auto recovered = state_from_rho_eps_ye(
        eos, changed.rho, changed.eps, changed.Ye);
    if (!std::isfinite(changed.kappa) ||
        fabs(changed.kappa - kappa) > 1.0e-10 * fmax(1.0, fabs(kappa)) ||
        fabs(recovered.temperature - state.temperature) >
            1.0e-7 * fmax(state.temperature, 1.0e-12))
      CCTK_ERROR("AsterSeeds test: density adjustment broke EOS closure");
  }
}
} // namespace

extern "C" void AsterSeeds_Test(CCTK_ARGUMENTS) {
  eos_3p_idealgas eos;
  eos_3p::range er{0.0, 1.0}, rr{1.0e-12, 1.0}, yr{0.0, 1.0};
  eos.init(2.0, 1.0, er, rr, yr);
  test_closure(&eos);
  if (global_eos_3p_tab3d)
    test_closure(global_eos_3p_tab3d);
  CCTK_INFO("AsterSeeds density-adjustment closure tests passed");
}
} // namespace AsterSeeds
