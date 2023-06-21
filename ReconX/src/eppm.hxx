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
#define MAX5(a, b, c, d, e) (MAX(a, MAX(b, MAX(c, MAX(d, e)))))

namespace ReconX {

using namespace std;
using namespace Arith;
using namespace Loop;

enum class interface_t { Minus = 0, Plus = 1 };

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
approx_at_cell_interface(const GF3D2<const CCTK_REAL> &gf,
                         const array<const vect<int, dim>, 5> &cells,
                         const CCTK_INT interface) {
  const auto &Im = cells.at(interface);
  const auto &I = cells.at(interface + 1);
  const auto &Ip = cells.at(interface + 2);
  const auto &Ipp = cells.at(interface + 3);

  return 7.0 / 12.0 * (gf(I) + gf(Ip)) - (1.0 / 12.0) * (gf(Im) + gf(Ipp));
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
limit(const GF3D2<const CCTK_REAL> &gf,
      const array<const vect<int, dim>, 5> &cells, const CCTK_REAL gf_rc,
      const CCTK_INT interface, const CCTK_REAL C) {
  const auto &Im = cells.at(interface);
  const auto &I = cells.at(interface + 1);
  const auto &Ip = cells.at(interface + 2);
  const auto &Ipp = cells.at(interface + 3);

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

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void
monotonize(const GF3D2<const CCTK_REAL> &gf,
           const array<const vect<int, dim>, 5> &cells, CCTK_REAL &rc_minus,
           CCTK_REAL &rc_plus, const CCTK_REAL C) {
  const auto &Imm = cells.at(0);
  const auto &Im = cells.at(1);
  const auto &I = cells.at(2);
  const auto &Ip = cells.at(3);
  const auto &Ipp = cells.at(4);

  CCTK_REAL D2aLim = 0;
  CCTK_REAL rhi = 0;
  const CCTK_REAL daplus = rc_plus - gf(I);
  const CCTK_REAL daminus = gf(I) - rc_minus;
  if (daplus * daminus <= 0 || (gf(Imm) - gf(I)) * (gf(I) - gf(Ipp)) <= 0) {
    const CCTK_REAL D2a = -(12.0 * gf(I) - 6.0 * (rc_minus + rc_plus));
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
        rc_plus = gf(I) + daplus * rhi;
        rc_minus = gf(I) - daminus * rhi;
      } else if (fabs(daminus) >= 2.0 * fabs(daplus)) {
        rc_minus = gf(I) - (2.0 * (1.0 - rhi) * daplus + rhi * daminus);
      } else if (fabs(daplus) >= 2.0 * fabs(daminus)) {
        rc_plus = gf(I) + (2.0 * (1.0 - rhi) * daminus + rhi * daplus);
      }
    }
  } else {
    if (fabs(daplus) >= 2.0 * fabs(daminus))
      rc_plus = gf(I) + 2.0 * daminus;
    else if (fabs(daminus) >= 2.0 * fabs(daplus))
      rc_minus = gf(I) - 2.0 * daplus;
  }

  return;
}

/* ePPM reconstruction scheme. (see Reisswig et al. 2013) based on McCorquodale
   & Colella (2011) */
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<CCTK_REAL, 2>
eppm(const GF3D2<const CCTK_REAL> &gf_var,
     const array<const vect<int, dim>, 5> &cells,
     const GF3D2<const CCTK_REAL> &gf_press,
     const GF3D2<const CCTK_REAL> &gf_vel_dir,
     const reconstruct_params_t &reconstruct_params) {
  // Unpack all cells in the stencil
  const auto &Imm = cells.at(0);
  const auto &Im = cells.at(1);
  const auto &I = cells.at(2);
  const auto &Ip = cells.at(3);
  const auto &Ipp = cells.at(4);

  // Unpack all PPM parameters
  const CCTK_REAL &ppm_eps = reconstruct_params.ppm_eps;
  const CCTK_REAL &ppm_small = reconstruct_params.ppm_small;
  const CCTK_REAL &ppm_omega1 = reconstruct_params.ppm_omega1;
  const CCTK_REAL &ppm_omega2 = reconstruct_params.ppm_omega2;
  const CCTK_REAL &enhanced_ppm_C2 = reconstruct_params.enhanced_ppm_C2;

  /* approx at cell interface */
  CCTK_REAL rc_minus =
      approx_at_cell_interface(gf_var, cells, CCTK_INT(interface_t::Minus));
  CCTK_REAL rc_plus =
      approx_at_cell_interface(gf_var, cells, CCTK_INT(interface_t::Plus));

  /* limit */
  rc_minus = limit(gf_var, cells, rc_minus, CCTK_INT(interface_t::Minus),
                   enhanced_ppm_C2);
  rc_plus = limit(gf_var, cells, rc_plus, CCTK_INT(interface_t::Plus),
                  enhanced_ppm_C2);

  /* monotonize */
  monotonize(gf_var, cells, rc_minus, rc_plus, enhanced_ppm_C2);

  /* apply flattening */
  const CCTK_REAL dpress_I = gf_press(Ip) - gf_press(Im);
  const CCTK_REAL dpress2 = gf_press(Ipp) - gf_press(Imm);
  const CCTK_REAL w_I =
      ((fabs(dpress_I) > ppm_eps * MIN(gf_press(Im), gf_press(Ip))) &&
       (gf_vel_dir(Im) > gf_vel_dir(Ip)))
          ? 1.0
          : 0.0;
  const CCTK_REAL ftilde_I =
      (fabs(dpress2) < ppm_small)
          ? 1.0
          : MAX(0.0, 1.0 - w_I * MAX(0.0, ppm_omega2 * (dpress_I / dpress2 -
                                                        ppm_omega1)));
  rc_minus = ftilde_I * rc_minus + (1.0 - ftilde_I) * gf_var(I);
  rc_plus = ftilde_I * rc_plus + (1.0 - ftilde_I) * gf_var(I);

  // Return the lower and upper reconstructed states in cell I
  return array<CCTK_REAL, 2>{rc_minus, rc_plus};
}

} // namespace ReconX

#endif // RECONX_EPPM_HXX
