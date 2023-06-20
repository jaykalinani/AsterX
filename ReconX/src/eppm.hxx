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
#define MIN4(a, b, c, d) (MIN(a, MIN(b, MIN(c, d))))

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
      const array<const vect<int, dim>, 5> &cells, const CCTK_REAL gf_rc,
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

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<CCTK_REAL, 2>
monotonize(const GF3D2<const CCTK_REAL> &gf,
           const array<const vect<int, dim>, 5> &cells,
           const array<CCTK_REAL, 2> gf_rc) {
  const auto &Imm = cells.at(0);
  const auto &Im = cells.at(1);
  const auto &I = cells.at(2);
  const auto &Ip = cells.at(3);
  const auto &Ipp = cells.at(4);

  CCTK_REAL D2aLim = 0;
  CCTK_REAL rhi = 0;
  const CCTK_REAL daplus = gf_rc[1] - gf(I);
  const CCTK_REAL daminus = gf(I) - gf_rc[0];
  if (daplus * daminus <= 0 || (gf(Imm) - gf(I)) * (gf(I) - gf(Ipp)) <= 0) {
    const CCTK_REAL D2a = -(12.0 * gf(I) - 6.0 * (gf_rc[0] + gf_rc[1]));
    const CCTK_REAL D2aC = (gf(Im) + gf(Ip)) - 2.0 * gf(I);
    const CCTK_REAL D2aL = (gf(Imm) + gf(I)) - 2.0 * gf(Im);
    const CCTK_REAL D2aR = (gf(I) + gf(Ipp)) - 2.0 * gf(Ip);
    if (copysign(1.0, D2a) == copysign(1.0, D2aC) &&
        copysign(1.0, D2a) == copysign(1.0, D2aL) &&
        copysign(1.0, D2a) == copysign(1.0, D2aR))
      D2aLim = copysign(1.0, D2a) *
               MIN4(C * fabs(D2aL), C * fabs(D2aR), C * fabs(D2aC), fabs(D2a));
    if (!(fabs(D2a) <= 1e-12 * MAX5(fabs(gf(Imm)), fabs(gf(Im)), fabs(gf(I)),
                                    fabs(gf(Ip)), fabs(gf(Ipp)))))
      rhi = D2aLim / D2a;
    if (!(rhi >= 1.0 - 1e-12)) {
      if (daplus * daminus < 0) {
        gf_rc_plus = gf(I) + daplus * rhi;
        gf_rc_minus = gf(I) - daminus * rhi;
      } else if (fabs(daminus) >= 2.0 * fabs(daplus)) {
        gf_rc_minus = gf(I) - (2.0 * (1.0 - rhi) * daplus + rhi * daminus);
      } else if (fabs(daplus) >= 2.0 * fabs(daminus)) {
        gf_rc_plus = gf(I) + (2.0 * (1.0 - rhi) * daminus + rhi * daplus);
      }
    }
  } else {
    if (fabs(daplus) >= 2.0 * fabs(daminus))
      gf_rc_plus = gf(I) + 2.0 * daminus;
    else if (fabs(daminus) >= 2.0 * fabs(daplus))
      gf_rc_minus = gf(I) - 2.0 * daplus;
  }

  return array<CCTK_REAL, 2>{gf_rc_minus, gf_rc_minus};
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
  {rc_minus, rc_plus} = monotonize(gf_var, cells, {rc_minus, rc_plus});

  /* apply flattening */
  const CCTK_REAL dpress_I = press_Ip - press_Im;
  const CCTK_REAL dpress2 = press_Ipp - press_Imm;
  const CCTK_REAL w_I = ((fabs(dpress_I) > ppm_eps * MIN(press_Im, press_Ip)) &&
                         (gf_vel_dir(Im) > gf_vel_dir(Ip)))
                            ? 1.0
                            : 0.0;
  const CCTK_REAL ftilde_I =
      (fabs(dpress2) < ppm_small)
          ? 1.0
          : MAX(0.0, 1.0. - w_I * MAX(0.0, ppm_omega2 * (dpress_I / dpress2 -
                                                         ppm_omega1)));
  const CCTK_REAL one_minus_ftilde_I_gfI = (1 - ftilde_I) * gf_I;
  rc_minus = ftilde_I * rc_minus + one_minus_ftilde_I_gfI;
  rc_plus = ftilde_I * rc_plus + one_minus_ftilde_I_gfI;

  // Return the lower and upper reconstructed states in cell I
  return array<CCTK_REAL, 2> {rc_minus, rc_plus};
}

} // namespace ReconX

#endif // RECONX_EPPM_HXX
