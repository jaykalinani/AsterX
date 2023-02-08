#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Arith;

inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST vec<vec<CCTK_REAL, 4>, 2>
    eigenvalues(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                vec<CCTK_REAL, 2> vel, vec<CCTK_REAL, 2> rho,
                vec<CCTK_REAL, 2> cs2, vec<CCTK_REAL, 2> w_lor,
                vec<CCTK_REAL, 2> h, vec<CCTK_REAL, 2> bsq) {
  // computing characteristics for the minus side
  // See Eq. (28) of Giacomazzo & Rezzolla (2007) with b^i=0
  vec<CCTK_REAL, 3> a_m{
      (bsq(0) + cs2(0) * h(0) * rho(0)) *
              (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
          (-1 + cs2(0)) * h(0) * rho(0) * pow2(beta_avg - alp_avg * vel(0)) *
              pow2(w_lor(0)),

      2 * beta_avg * (bsq(0) + cs2(0) * h(0) * rho(0)) -
          2 * (-1 + cs2(0)) * h(0) * rho(0) * (beta_avg - alp_avg * vel(0)) *
              pow2(w_lor(0)),

      bsq(0) +
          h(0) * rho(0) * (cs2(0) + pow2(w_lor(0)) - cs2(0) * pow2(w_lor(0)))};

  CCTK_REAL det_m = pow2(a_m(1)) - 4 * a_m(2) * a_m(0);
  if (det_m < 0)
    det_m = 0;

  vec<CCTK_REAL, 4> lambda_m{((-a_m(1) + sqrt(det_m)) / (2 * a_m(2))) / alp_avg,
                             ((-a_m(1) + sqrt(det_m)) / (2 * a_m(2))) / alp_avg,
                             ((-a_m(1) - sqrt(det_m)) / (2 * a_m(2))) / alp_avg,
                             ((-a_m(1) - sqrt(det_m)) / (2 * a_m(2))) /
                                 alp_avg};

  // computing characteristics for the plus side

  vec<CCTK_REAL, 3> a_p{
      (bsq(1) + cs2(1) * h(1) * rho(1)) *
              (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
          (-1 + cs2(1)) * h(1) * rho(1) * pow2(beta_avg - alp_avg * vel(1)) *
              pow2(w_lor(1)),

      2 * beta_avg * (bsq(1) + cs2(1) * h(1) * rho(1)) -
          2 * (-1 + cs2(1)) * h(1) * rho(1) * (beta_avg - alp_avg * vel(1)) *
              pow2(w_lor(1)),

      bsq(1) +
          h(1) * rho(1) * (cs2(1) + pow2(w_lor(1)) - cs2(1) * pow2(w_lor(1)))};

  CCTK_REAL det_p = pow2(a_p(1)) - 4 * a_p(2) * a_p(0);
  if (det_p < 0)
    det_p = 0;

  vec<CCTK_REAL, 4> lambda_p{((-a_p(1) + sqrt(det_p)) / (2 * a_p(2))) / alp_avg,
                             ((-a_p(1) + sqrt(det_p)) / (2 * a_p(2))) / alp_avg,
                             ((-a_p(1) - sqrt(det_p)) / (2 * a_p(2))) / alp_avg,
                             ((-a_p(1) - sqrt(det_p)) / (2 * a_p(2))) /
                                 alp_avg};

  // 2D array containing characteristics for left (minus) and right
  // (plus) sides
  vec<vec<CCTK_REAL, 4>, 2> lambda{lambda_m, lambda_p};
  return lambda;
};

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
laxf(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
     vec<CCTK_REAL, 2> flux) {
  const CCTK_REAL charmax =
      max({CCTK_REAL(0), fabs(lam(0)(0)), fabs(lam(0)(1)), fabs(lam(0)(2)),
           fabs(lam(0)(3)), fabs(lam(1)(0)), fabs(lam(1)(1)), fabs(lam(1)(2)),
           fabs(lam(1)(3))});

  return 0.5 * ((flux(0) + flux(1)) - charmax * (var(1) - var(0)));
}

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST CCTK_REAL
hlle(vec<vec<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
     vec<CCTK_REAL, 2> flux) {
  const CCTK_REAL charmax =
      max({CCTK_REAL(0), lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
           lam(1)(1), lam(1)(2), lam(1)(3)});

  const CCTK_REAL charmin =
      min({CCTK_REAL(0), lam(0)(0), lam(0)(1), lam(0)(2), lam(0)(3), lam(1)(0),
           lam(1)(1), lam(1)(2), lam(1)(3)});

  const CCTK_REAL charpm = charmax - charmin;

  return (charmax * flux(1) - charmin * flux(0) +
          charmax * charmin * (var(1) - var(0))) /
         charpm;
}

} // namespace AsterX
