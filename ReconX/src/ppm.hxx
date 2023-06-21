#ifndef RECONX_PPM_HXX
#define RECONX_PPM_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>

#include "reconx_utils.hxx"

namespace ReconX {

using namespace std;
using namespace Arith;
using namespace Loop;

/* PPM reconstruction scheme. See Colella & Woodward (1984) (e.g. at
 * https://crd.lbl.gov/assets/pubs_presos/AMCS/ANAG/A141984.pdf) */
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<CCTK_REAL, 2>
ppm(const GF3D2<const CCTK_REAL> &gf_var,
    const array<const vect<int, dim>, 5> &cells, const CCTK_INT &dir,
    const bool &gf_is_rho, const GF3D2<const CCTK_REAL> &gf_press,
    const GF3D2<const CCTK_REAL> &gf_vel_dir,
    const reconstruct_params_t &reconstruct_params) {
  // Unpack all cells in the stencil
  const auto &Imm = cells[0];
  const auto &Im = cells[1];
  const auto &I = cells[2];
  const auto &Ip = cells[3];
  const auto &Ipp = cells[4];

  // Unpack all PPM parameters
  const bool &ppm_shock_detection = reconstruct_params.ppm_shock_detection;
  const bool &ppm_zone_flattening = reconstruct_params.ppm_zone_flattening;
  const CCTK_REAL &poly_k = reconstruct_params.poly_k;
  const CCTK_REAL &poly_gamma = reconstruct_params.poly_gamma;
  const CCTK_REAL &ppm_eta1 = reconstruct_params.ppm_eta1;
  const CCTK_REAL &ppm_eta2 = reconstruct_params.ppm_eta2;
  const CCTK_REAL &ppm_eps = reconstruct_params.ppm_eps;
  const CCTK_REAL &ppm_eps_shock = reconstruct_params.ppm_eps_shock;
  const CCTK_REAL &ppm_small = reconstruct_params.ppm_small;
  const CCTK_REAL &ppm_omega1 = reconstruct_params.ppm_omega1;
  const CCTK_REAL &ppm_omega2 = reconstruct_params.ppm_omega2;

  // Grid function at neighboring cells
  const CCTK_REAL &gf_Imm = gf_var(Imm);
  const CCTK_REAL &gf_Im = gf_var(Im);
  const CCTK_REAL &gf_I = gf_var(I);
  const CCTK_REAL &gf_Ip = gf_var(Ip);
  const CCTK_REAL &gf_Ipp = gf_var(Ipp);

  // Helpers
  const CCTK_REAL diff_Im = gf_I - gf_Imm;
  const CCTK_REAL diff_I = gf_Ip - gf_Im;
  const CCTK_REAL diff_Ip = gf_Ipp - gf_I;

  const CCTK_REAL delta_Im = 0.5 * diff_Im;
  const CCTK_REAL delta_I = 0.5 * diff_I;
  const CCTK_REAL delta_Ip = 0.5 * diff_Ip;

  const CCTK_REAL _2fabs_Im_Imm = 2 * fabs(gf_Im - gf_Imm);
  const CCTK_REAL _2fabs_I_Im = 2 * fabs(gf_I - gf_Im);
  const CCTK_REAL _2fabs_Ip_I = 2 * fabs(gf_Ip - gf_I);
  const CCTK_REAL _2fabs_Ipp_Ip = 2 * fabs(gf_Ipp - gf_Ip);

  const bool same_sgn_Im = ((gf_I - gf_Im) * (gf_Im - gf_Imm) > 0);
  const bool same_sgn_I = ((gf_Ip - gf_I) * (gf_I - gf_Im) > 0);
  const bool same_sgn_Ip = ((gf_Ipp - gf_Ip) * (gf_Ip - gf_I) > 0);

  const CCTK_REAL deltamod_Im =
      same_sgn_Im
          ? sgn(delta_Im) * min(fabs(delta_Im), min(_2fabs_Im_Imm, _2fabs_I_Im))
          : 0;
  const CCTK_REAL deltamod_I =
      same_sgn_I
          ? sgn(delta_I) * min(fabs(delta_I), min(_2fabs_I_Im, _2fabs_Ip_I))
          : 0;
  const CCTK_REAL deltamod_Ip =
      same_sgn_Ip
          ? sgn(delta_Ip) * min(fabs(delta_Ip), min(_2fabs_Ip_I, _2fabs_Ipp_Ip))
          : 0;

  /* Initial reconstructed states at the interfaces between cells Im/I ans I/Ip
   * NOTE: not const because they may change later */
  const CCTK_REAL gf_Im_I =
      0.5 * (gf_Im + gf_I) + (1. / 6.) * (deltamod_Im - deltamod_I);
  const CCTK_REAL gf_I_Ip =
      0.5 * (gf_I + gf_Ip) + (1. / 6.) * (deltamod_I - deltamod_Ip);

  CCTK_REAL rc_low = gf_Im_I;
  CCTK_REAL rc_up = gf_I_Ip;

  /* Pressure may or may not be needed, or may not be needed at all points
   * (depending on whether shock detection and zone flattening are activated or
   * not) */
  // const CCTK_REAL &press_Immm = gf_press(Immm);  // Only used in the original
  // zone flattening scheme
  const CCTK_REAL &press_Imm = gf_press(Imm);
  const CCTK_REAL &press_Im = gf_press(Im);
  // const CCTK_REAL &press_I    = gf_press(I);     // Only used in the original
  // zone flattening scheme
  const CCTK_REAL &press_Ip = gf_press(Ip);
  const CCTK_REAL &press_Ipp = gf_press(Ipp);
  // const CCTK_REAL &press_Ippp = gf_press(Ippp);  // Only used in the original
  // zone flattening scheme

  const CCTK_REAL diff_press_I = press_Ip - press_Im;
  const CCTK_REAL min_press_I = min(press_Im, press_Ip);

  /* Shock detection (eqs. 1.15, 1.16, 1.17 with uniform grid spacing).
   * This is only applied to rho and only if the shock is marked as a contact
   * discontinuity (see eq. 3.2) */
  // FIXME: contact discontinuity check only valid for polytropic/ideal-fluid
  // EOS
  if (ppm_shock_detection and gf_is_rho) {
    const CCTK_REAL k_gamma = poly_k * poly_gamma;
    const bool contact_disc = (k_gamma * fabs(diff_I) * min_press_I >=
                               fabs(diff_press_I) * min(gf_Ip, gf_Im));
    if (contact_disc) {
      // This assumes gf_var is rho (gf_is_rho is true)
      const CCTK_REAL d2rho_Im = gf_I - 2 * gf_Im + gf_Imm;
      const CCTK_REAL d2rho_Ip = gf_Ipp - 2 * gf_Ip + gf_I;
      const CCTK_REAL d2_prod = d2rho_Im * d2rho_Ip;
      const bool cond2 = (fabs(diff_I) - ppm_eps_shock*min(fabs(gf_Ip), fabs(gf_Im))) > 0.;
      
      const CCTK_REAL eta_tilde_I =
          ((d2_prod < 0) and cond2) ? ((-1. / 6.) * (d2rho_Ip - d2rho_Im) / diff_I) : 0;
      const CCTK_REAL eta_I =
          max(0., min(ppm_eta1 * (eta_tilde_I - ppm_eta2), 1.));

      rc_low = (1 - eta_I) * rc_low + eta_I * (gf_Im + 0.5 * deltamod_Im);
      rc_up = (1 - eta_I) * rc_up + eta_I * (gf_Ip - 0.5 * deltamod_Ip);
    }
  }

  /* Zone flattening to reduce post-shock oscillations (see eq. 4.1 + appendix).
   * The flattening parameter f_I is set to ftilde_I instead of
   * max(ftilde_I, ftilde_{I+s_I}) (where s_I can be +1 or -1), thereby avoiding
   * using four ghost cells at interprocess/domain boundaries. This should not
   * be a major issue and is done in many GRMHD codes (e.g. WhiskyMHD, GRHydro,
   * Spritz, IllinoisGRMHD). */
  if (ppm_zone_flattening) {
    const CCTK_REAL w_I = ( (fabs(diff_press_I) > ppm_eps * min_press_I) and
                           (gf_vel_dir(Im) > gf_vel_dir(Ip)))
                              ? 1.
                              : 0;
    const CCTK_REAL ftilde_I = (fabs(press_Ipp - press_Imm) < ppm_small) ? 1.0 :
         max(0., 1. - w_I * max(0., ppm_omega2 *
                        ( (diff_press_I / (press_Ipp - press_Imm)) - ppm_omega1)));

    const CCTK_REAL one_minus_ftilde_I_gfI = (1 - ftilde_I)*gf_I;

    rc_low = ftilde_I * rc_low + one_minus_ftilde_I_gfI;
    rc_up = ftilde_I * rc_up + one_minus_ftilde_I_gfI;

    // This would require one more ghost cell and it's not worth it
    /*if (diff_press_I < 0) {
      const CCTK_REAL diff_press_Ip = press_Ipp - press_I;
      const CCTK_REAL w_Ip          = (diff_press_Ip > ppm_eps*min(press_I,
    press_Ipp) and gf_vel_dir(I) > gf_vel_dir(Ipp)) ? 1 : 0; const CCTK_REAL
    ftilde_Ip     = 1 - max(0., w_Ip*ppm_omega2*(diff_press_Ip/(press_Ippp -
    press_Im) - ppm_omega1)); const CCTK_REAL f_I           = max(ftilde_I,
    ftilde_Ip); const CCTK_REAL one_minus_fI  = 1 - f_I; const CCTK_REAL fI_gfI
    = f_I*gf_I; rc_low = fI_gfI + one_minus_fI*rc_low; rc_up  = fI_gfI +
    one_minus_fI*rc_up;
    }

    else {
      const CCTK_REAL diff_press_Im = press_I - press_Imm;
      const CCTK_REAL w_Im          = (diff_press_Im > ppm_eps*min(press_Imm,
    press_I) and gf_vel_dir(Imm) > gf_vel_dir(I)) ? 1 : 0; const CCTK_REAL
    ftilde_Im     = 1 - max(0., w_Im*ppm_omega2*(diff_press_Im/(press_Ip -
    press_Immm) - ppm_omega1)); const CCTK_REAL f_I           = max(ftilde_I,
    ftilde_Im); const CCTK_REAL one_minus_fI  = 1 - f_I; const CCTK_REAL fI_gfI
    = f_I*gf_I; rc_low = fI_gfI + one_minus_fI*rc_low; rc_up  = fI_gfI +
    one_minus_fI*rc_up;
    }*/
  }

  // Monotonization (see eq. 1.10)
  if ((rc_up - gf_I) * (gf_I - rc_low) <= 0) {
    rc_low = gf_I;
    rc_up = gf_I;
  }
  else {
    const CCTK_REAL diff_rc = rc_up - rc_low;
    const CCTK_REAL diff_rc_sq = diff_rc * diff_rc;
    const CCTK_REAL gf6_I = 6 * diff_rc * (gf_I - 0.5 * (rc_low + rc_up));

    if (gf6_I > diff_rc_sq) {
      rc_low = 3 * gf_I - 2 * rc_up;
    }

    else if (gf6_I < -diff_rc_sq) {
      rc_up = 3 * gf_I - 2 * rc_low;
    }
  }

  // Return the lower and upper reconstructed states in cell I
  const array<CCTK_REAL, 2> rc = {rc_low, rc_up};
  return rc;

}

} // namespace ReconX

#endif // RECONX_PPM_HXX
