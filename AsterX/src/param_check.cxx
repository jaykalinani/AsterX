#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cmath>

#include "atmo.hxx"
#include "atmo_cache.hxx"
#include "setup_eos.hxx"

namespace AsterX {
using namespace EOSX;
using namespace Con2PrimFactory;

template <typename EOSIDType, typename EOSType>
void ReportAtmosphere(const EOSIDType *eos_1p, const EOSType *eos_3p) {
  DECLARE_CCTK_PARAMETERS;
  const auto atmo = make_atmosphere(
      eos_1p, eos_3p, 0.0, rho_abs_min, p_atmo, t_atmo, Ye_atmo, r_atmo,
      n_rho_atmo, n_press_atmo, n_temp_atmo, atmo_tol, thermal_eos_atmo,
      use_press_atmo);
  if (!std::isfinite(atmo.rho_atmo) || !std::isfinite(atmo.eps_atmo) ||
      !std::isfinite(atmo.press_atmo) || !std::isfinite(atmo.temp_atmo) ||
      !std::isfinite(atmo.entropy_atmo) || atmo.rho_atmo <= 0.0)
    CCTK_ERROR("The selected atmosphere does not produce a finite EOS state");
  CCTK_VINFO("Inner atmosphere: rho=%.16e eps=%.16e P=%.16e T=%.16e "
             "Ye=%.16e kappa=%.16e", atmo.rho_atmo, atmo.eps_atmo,
             atmo.press_atmo, atmo.temp_atmo, atmo.ye_atmo,
             atmo.entropy_atmo);
  if (atmo.rho_atmo != rho_abs_min || atmo.ye_atmo != Ye_atmo)
    CCTK_INFO("Requested atmosphere rho or Ye was limited to the EOS domain");
  cached_atmo.store(eos_1p, eos_3p, atmo, rho_abs_min, p_atmo, t_atmo,
                     Ye_atmo, thermal_eos_atmo, use_press_atmo);
}

extern "C" void AsterX_ParamCheck(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_ParamCheck;
  DECLARE_CCTK_PARAMETERS;

  // CarpetX runs PARAMCHECK after EOS setup and before initial data,
  // including on restart. Check parameter combinations only here.
  cached_atmo.valid = false;
  if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    if (!use_temperature || !reconstruct_with_temperature ||
        !thermal_eos_atmo || use_press_atmo)
      CCTK_ERROR("Tabulated3d requires use_temperature=yes, "
                 "reconstruct_with_temperature=yes, thermal_eos_atmo=yes "
                 "and use_press_atmo=no");
    if (use_entropy_fix || CCTK_EQUALS(c2p_prime, "Entropy") ||
        CCTK_EQUALS(c2p_second, "Entropy"))
      CCTK_ERROR("Tabulated3d does not implement the inversions needed by "
                 "the entropy C2P; disable use_entropy_fix and Entropy solvers");
  }
  if (CCTK_EQUALS(evolution_eos, "Hybrid") && thermal_eos_atmo)
    CCTK_ERROR("Hybrid EOS requires thermal_eos_atmo=no");
  if (!(sigma_max > 0.0) || !(inv_beta_max > 0.0) || !(c2p_tol > 0.0))
    CCTK_ERROR("sigma_max, inv_beta_max and c2p_tol must be positive");

  if (local_spatial_order != 2 && local_spatial_order != 4)
    CCTK_ERROR("local_spatial_order must be set to 2 or 4");
  if (tmunu_interp_order != 2 && tmunu_interp_order != 4)
    CCTK_ERROR("tmunu_interp_order must be set to 2 or 4");

  if (!freeze_evolution) {
    const auto nghost = [](const char *method) {
      if (CCTK_EQUALS(method, "Godunov"))
        return 1;
      if (CCTK_EQUALS(method, "minmod") ||
          CCTK_EQUALS(method, "monocentral"))
        return 2;
      return 3;
    };
    int need = std::max(nghost(reconstruction_method), nghost(loworder_method));
    if (add_dissipation)
      need = std::max(need, 3);
    for (int d = 0; d < 3; ++d)
      if (cctk_nghostzones[d] < need)
        CCTK_VERROR("Reconstruction and dissipation require >=%d ghost zones; "
                   "have (%d,%d,%d).", need, cctk_nghostzones[0],
                   cctk_nghostzones[1], cctk_nghostzones[2]);
  }

  if (!thermal_eos_atmo) {
    CCTK_INFO("Atmosphere mode: cold initial-data EOS energy, closed with "
              "the evolution EOS. t_atmo, p_atmo and their grading "
              "exponents do not construct the atmosphere in this mode");
    if (CCTK_EQUALS(evolution_eos, "IdealGas") &&
        CCTK_EQUALS(initial_data_eos, "Polytropic") && poly_gamma != gl_gamma)
      CCTK_WARN(CCTK_WARN_ALERT,
                "poly_gamma differs from gl_gamma: cold energy is matched, "
                "but cold and evolution pressures will differ");
  } else if (use_press_atmo) {
    CCTK_INFO("Atmosphere mode: pressure-primary; p_atmo and n_press_atmo "
              "are active, t_atmo and n_temp_atmo are ignored");
  } else {
    CCTK_INFO("Atmosphere mode: temperature-primary; t_atmo and n_temp_atmo "
              "are active, p_atmo and n_press_atmo are ignored");
  }
  CCTK_INFO("eps_atmo is a legacy input; atmosphere energy is derived "
            "from the selected EOS state");
  CCTK_VINFO("Atmosphere grading: r_atmo=%.16e n_rho=%.16e n_temp=%.16e "
             "n_press=%.16e", r_atmo, n_rho_atmo, n_temp_atmo, n_press_atmo);
  const CCTK_REAL face_atmo_factor =
      recon_use_atmo_tol ? 1.0 + atmo_tol : recon_thresh;
  CCTK_VINFO("Atmosphere cutoff factors: cells=%.16e faces=%.16e; each "
             "multiplies its own local graded atmosphere density",
             1.0 + atmo_tol, face_atmo_factor);
  CCTK_VINFO("Conservative tau repair margin=%.16e; applied only after "
             "the global physical energy admissibility test fails",
             tauFluid_atmo);

  // Dispatch the same cold EOS that is used by the atmosphere builder.
  const auto report = [&](const auto *eos_3p) {
    if (global_eos_1p_pwpoly)
      ReportAtmosphere(global_eos_1p_pwpoly, eos_3p);
    else if (global_eos_1p_poly)
      ReportAtmosphere(global_eos_1p_poly, eos_3p);
    else
      CCTK_ERROR("No initial-data EOS was initialized for the atmosphere");
  };
  if (global_eos_3p_ig)
    report(global_eos_3p_ig);
  else if (global_eos_3p_tab3d)
    report(global_eos_3p_tab3d);
  // Hybrid is outside this repair series; its existing setup remains active.
}

} // namespace AsterX
