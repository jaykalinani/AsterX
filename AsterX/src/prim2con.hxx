#ifndef PRIM2CON_HXX
#define PRIM2CON_HXX

#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>
#include "utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

struct prim {
  CCTK_REAL rho;
  CCTK_REAL velx, vely, velz;
  CCTK_REAL eps, press;
  CCTK_REAL Bvecx, Bvecy, Bvecz;
};

struct cons {
  CCTK_REAL dens;
  CCTK_REAL momx, momy, momz;
  CCTK_REAL tau;
  CCTK_REAL dBvecx, dBvecy, dBvecz;
};

CCTK_DEVICE CCTK_HOST void prim2con(const smat<CCTK_REAL, 3> &g,
                                    const CCTK_REAL &lapse,
                                    const vec<CCTK_REAL, 3> &beta_up,
                                    const prim &pv, cons &cv) {

  // determinant of spatial metric
  const CCTK_REAL sqrt_detg = sqrt(calc_det(g));

  // TODO: compute specific internal energy based on user-specified EOS
  // currently, computing eps for classical ideal gas

  /* Computing v_j */
  const vec<CCTK_REAL, 3> v_up{pv.velx, pv.vely, pv.velz};
  const vec<CCTK_REAL, 3> v_low = calc_contraction(g, v_up);

  const CCTK_REAL w_lorentz = calc_wlorentz(v_low, v_up);

  /* Computing beta_j */
  const vec<CCTK_REAL, 3> beta_low = calc_contraction(g, beta_up);

  /* Computing B_j */
  const vec<CCTK_REAL, 3> B_up{pv.Bvecx, pv.Bvecy, pv.Bvecz};
  const vec<CCTK_REAL, 3> B_low = calc_contraction(g, B_up);

  /* Computing b^t : this is b^0 * alp */
  const CCTK_REAL bst = w_lorentz * calc_contraction(B_up, v_low);

  /* Computing b^j */
  const vec<CCTK_REAL, 3> b_up([&](int i) ARITH_INLINE {
    return (B_up(i) / w_lorentz + bst * (v_up(i) - beta_up(i) / lapse));
  });

  /* Computing b_j */
  const vec<CCTK_REAL, 3> gb = calc_contraction(g, b_up);
  const vec<CCTK_REAL, 3> b_low([&](int i) ARITH_INLINE {
    return gb(i) + beta_low(i) * bst / lapse;
  });

  /* Computing b^mu b_mu */
  const CCTK_REAL bs2 =
      (calc_contraction(B_up, B_low) + bst * bst) / (w_lorentz * w_lorentz);

  // computing conserved from primitives
  cv.dens = sqrt_detg * pv.rho * w_lorentz;

  cv.momx =
      sqrt_detg * (w_lorentz * w_lorentz *
                       (pv.rho * (1.0 + pv.eps) + pv.press + bs2) * v_low(0) -
                   bst * b_low(0));

  cv.momy =
      sqrt_detg * (w_lorentz * w_lorentz *
                       (pv.rho * (1.0 + pv.eps) + pv.press + bs2) * v_low(1) -
                   bst * b_low(1));

  cv.momz =
      sqrt_detg * (w_lorentz * w_lorentz *
                       (pv.rho * (1.0 + pv.eps) + pv.press + bs2) * v_low(2) -
                   bst * b_low(2));

  cv.tau = sqrt_detg * (w_lorentz * w_lorentz *
                            (pv.rho * (1.0 + pv.eps) + pv.press + bs2) -
                        (pv.press + 0.5 * bs2) - bst * bst) -
           cv.dens;

  cv.dBvecx = sqrt_detg * pv.Bvecx;
  cv.dBvecy = sqrt_detg * pv.Bvecy;
  cv.dBvecz = sqrt_detg * pv.Bvecz;
}

} // namespace AsterX

#endif // #ifndef PRIM2CON_HXX
