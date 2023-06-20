#ifndef RECONX_EPPM_HXX
#define RECONX_EPPM_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>

#include "reconx_utils.hxx"

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN3(a, b, c) (MIN(a, MIN(b, c)))

namespace ReconX {

using namespace std;
using namespace Arith;
using namespace Loop;

enum class interface_t { Minus = 0, Plus = 1 };

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
approx_at_cell_interface(const GF3D2<const CCTK_REAL> &gf,
                         const array<const vect<int, dim>, 5> &cells,
                         interface_t &interface) {
  const auto &Im = cells.at(CCTK_INT(interface));
  const auto &I = cells.at(CCTK_INT(interface) + 1);
  const auto &Ip = cells.at(CCTK_INT(interface) + 2);
  const auto &Ipp = cells.at(CCTK_INT(interface) + 3);

  return 7.0 / 12.0 * (gf(I) + gf(Ip)) - (1.0 / 12.0) * (gf(Im) + gf(Ipp));
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
limit(const GF3D2<const CCTK_REAL> &gf,
      const array<const vect<int, dim>, 5> &cells, CCTK_REAL gf_rc,
      interface_t &interface) {
  const auto &Im = cells.at(CCTK_INT(interface));
  const auto &I = cells.at(CCTK_INT(interface) + 1);
  const auto &Ip = cells.at(CCTK_INT(interface) + 2);
  const auto &Ipp = cells.at(CCTK_INT(interface) + 3);

  if ((MIN(gf(I), gf(Ip)) <= gf_rc) && (gf_rc <= MAX(gf(I), gf(Ip)))) {
    return gf_rc;
  } else {
    const CCTK_REAL D2a = 3.0 * ((gf(I) + gf(Ip)) - 2.0 * gf_rc);
    const CCTK_REAL D2aL = ((gf(Im) + gf(Ip)) - 2.0 * gf(I));
    const CCTK_REAL D2aR = ((gf(I) + gf(Ipp)) - 2.0 * gf(Ip));
    const CCTK_REAL D2aLim = copysign(1.0, D2a) *
                             MIN3(C * fabs(D2aL), C * fabs(D2aR), fabs(D2a)) *
                             1.0 / 3.0;
    if (D2a * D2aR >= 0 && D2a * D2aL >= 0)
      return 0.5 * (gf(I) + gf(Ip)) - D2aLim;
    else
      return 0.5 * (gf(I) + gf(Ip));
  }
}

/* ePPM reconstruction scheme. (see Reisswig et al. 2013) based on McCorquodale
   & Colella (2011) */
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<CCTK_REAL, 2>
eppm(const GF3D2<const CCTK_REAL> &gf_var,
     const array<const vect<int, dim>, 5> &cells, const CCTK_INT &dir,
     const bool &gf_is_rho, const GF3D2<const CCTK_REAL> &gf_press,
     const GF3D2<const CCTK_REAL> &gf_vel_dir,
     const reconstruct_params_t &reconstruct_params) {
  // Unpack all cells in the stencil
  const auto &Imm = cells.at(0);
  const auto &Im = cells.at(1);
  const auto &I = cells.at(2);
  const auto &Ip = cells.at(3);
  const auto &Ipp = cells.at(4);

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

  /* approx at cell interface */
  CCTK_REAL rc_minus =
      approx_at_cell_interface(gf_var, cells, interface_t::Minus);
  CCTK_REAL rc_plus =
      approx_at_cell_interface(gf_var, cells, interface_t::Plus);

  /* limit */
  rc_minus = limit(gf_var, cells, rc_minus, interface_t::Minus);
  rc_plus = limit(gf_var, cells, rc_minus, interface_t::Plus);

  /* monotonize */

  /* apply flattening */

  // Return the lower and upper reconstructed states in cell I
  const array<CCTK_REAL, 2> rc = {rc_minus, rc_plus};
  return rc;
}

} // namespace ReconX

#endif // RECONX_EPPM_HXX
