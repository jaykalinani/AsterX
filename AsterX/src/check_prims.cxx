#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include "aster_utils.hxx"
#include "atmo.hxx"
#include "setup_eos.hxx"

namespace AsterX {
using namespace AsterUtils;
using namespace Loop;
using namespace EOSX;
using namespace Con2PrimFactory;
using namespace std;

enum class eos_3param { IdealGas, Hybrid, Tabulated };

template <typename EOSIDType, typename EOSType>
void CheckPrims(CCTK_ARGUMENTS, EOSIDType *eos_1p, EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CheckPrims;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  grid.loop_all_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                   calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                   calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};

        vec<CCTK_REAL, 3> v_up{velx(p.I), vely(p.I), velz(p.I)};

        CCTK_REAL rhoL = rho(p.I);
        CCTK_REAL epsL = eps(p.I);
        CCTK_REAL pressL = press(p.I);
        CCTK_REAL YeL = Ye(p.I);
        CCTK_REAL tempL = temperature(p.I);
        CCTK_REAL entropyL;

        const CCTK_REAL radial_distance =
            sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
        const auto atmo = make_atmosphere(
            eos_1p, eos_3p, radial_distance, rho_abs_min, p_atmo, t_atmo,
            Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo, n_temp_atmo,
            atmo_tol, thermal_eos_atmo, use_press_atmo);

        const CCTK_REAL w_lim = sqrt(1.0 + vw_lim * vw_lim);
        const CCTK_REAL v_lim = vw_lim / w_lim;

        // ----------
        // Ceiling for velocity
        // ----------

        // check if computed velocities are within the specified limit
        vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);
        CCTK_REAL vsq_Sol = calc_contraction(v_low, v_up);
        CCTK_REAL sol_v = sqrt(vsq_Sol);
        if (sol_v > v_lim) {

          v_up *= v_lim / sol_v;
        }

        const bool finite = std::isfinite(rhoL) && std::isfinite(YeL) &&
            std::isfinite(use_temperature ? tempL :
                          (use_press_atmo ? pressL : epsL)) &&
            std::isfinite(v_up(0)) && std::isfinite(v_up(1)) &&
            std::isfinite(v_up(2));
        if (!finite || rhoL <= atmo.rho_cut) {
          // Reset the complete primitive state instead of retaining thermal
          // quantities from a cell that has been classified as atmosphere.
          rhoL = atmo.rho_atmo;
          epsL = atmo.eps_atmo;
          pressL = atmo.press_atmo;
          YeL = atmo.ye_atmo;
          tempL = atmo.temp_atmo;
          entropyL = atmo.entropy_atmo;
          v_up(0) = 0.0;
          v_up(1) = 0.0;
          v_up(2) = 0.0;
        } else {
          EOSX::thermo_state state;
          if (use_temperature) {
            tempL = fmax(tempL, atmo.temp_atmo);
            state = EOSX::state_from_rho_temp_ye(eos_3p, rhoL, tempL, YeL);
          } else if (use_press_atmo) {
            pressL = fmax(pressL, atmo.press_atmo);
            state = EOSX::state_from_rho_press_ye(eos_3p, rhoL, pressL, YeL);
          } else {
            state = EOSX::state_from_rho_eps_ye(eos_3p, rhoL, epsL, YeL);
            if (state.temperature < atmo.temp_atmo)
              state = EOSX::state_from_rho_temp_ye(
                  eos_3p, state.rho, atmo.temp_atmo, state.Ye);
          }

          rhoL = state.rho;
          epsL = state.eps;
          pressL = state.press;
          YeL = state.Ye;
          tempL = state.temperature;
          entropyL = state.kappa;
        }

        // ---------- End of validity check

        rho(p.I) = rhoL;
        velx(p.I) = v_up(0);
        vely(p.I) = v_up(1);
        velz(p.I) = v_up(2);
        eps(p.I) = epsL;
        press(p.I) = pressL;
        entropy(p.I) = entropyL;
        Ye(p.I) = YeL;
        temperature(p.I) = tempL;

        saved_rho(p.I) = rhoL;
        saved_velx(p.I) = v_up(0);
        saved_vely(p.I) = v_up(1);
        saved_velz(p.I) = v_up(2);
        saved_eps(p.I) = epsL;
        saved_Ye(p.I) = YeL;

        v_low = calc_contraction(g, v_up);
        CCTK_REAL wlor = calc_wlorentz(v_low, v_up);

        zvec_x(p.I) = wlor * v_up(0);
        zvec_y(p.I) = wlor * v_up(1);
        zvec_z(p.I) = wlor * v_up(2);

        svec_x(p.I) = (rhoL + rhoL * epsL + pressL) * wlor * wlor * v_up(0);
        svec_y(p.I) = (rhoL + rhoL * epsL + pressL) * wlor * wlor * v_up(1);
        svec_z(p.I) = (rhoL + rhoL * epsL + pressL) * wlor * wlor * v_up(2);
      });
}

extern "C" void AsterX_CheckPrims(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_CheckPrims;
  DECLARE_CCTK_PARAMETERS;

  eos_3param eos_3p_type;

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_3p_type = eos_3param::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_3p_type = eos_3param::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_3p_type = eos_3param::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_3p_type) {
  case eos_3param::IdealGas: {
    // Get local eos object
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_ig = global_eos_3p_ig;

    if (global_eos_1p_pwpoly)
      CheckPrims(cctkGH, global_eos_1p_pwpoly, eos_3p_ig);
    else
      CheckPrims(cctkGH, eos_1p_poly, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    if (global_eos_3p_hyb_pwpoly) {
      // pwpoly cold + pwpoly-hybrid
      auto eos_cold = global_eos_1p_pwpoly;
      auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;

      if (!eos_cold) {
        CCTK_ERROR("Hybrid(PWPolytropic) selected but no pwpoly cold EOS was "
                   "initialized");
      }
      CheckPrims(cctkGH, eos_cold, eos_3p_hyb);

    } else if (global_eos_3p_hyb_poly) {
      // poly cold + poly-hybrid
      auto eos_cold = global_eos_1p_poly;
      auto eos_3p_hyb = global_eos_3p_hyb_poly;

      if (!eos_cold) {
        CCTK_ERROR("Hybrid(Polytropic) selected but no polytropic cold EOS was "
                   "initialized");
      }
      CheckPrims(cctkGH, eos_cold, eos_3p_hyb);

    } else {
      CCTK_ERROR(
          "Hybrid EOS selected but no hybrid EOS object was initialized");
    }

    break;
  }
  case eos_3param::Tabulated: {
    // Get local eos object
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_tab3d = global_eos_3p_tab3d;

    if (global_eos_1p_pwpoly)
      CheckPrims(cctkGH, global_eos_1p_pwpoly, eos_3p_tab3d);
    else
      CheckPrims(cctkGH, eos_1p_poly, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}

} // namespace AsterX
