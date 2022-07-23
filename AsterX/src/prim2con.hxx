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

struct metric {
  CCTK_REAL gxx, gxy, gxz, gyy, gyz, gzz;
};

struct shift {
  CCTK_REAL betax, betay, betaz;
};

struct lapse {
  CCTK_REAL alp;
};

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

CCTK_DEVICE CCTK_HOST void prim2con(const metric &g, const lapse &lap,
                                    const shift &shft, const prim &pv,
                                    cons &cv) {

  // determinant of spatial metric
  const CCTK_REAL detg = calc_detg(g.gxx, g.gxy, g.gxz, g.gyy, g.gyz, g.gzz);
  const CCTK_REAL sqrt_detg = sqrt(detg);

  // TODO: compute specific internal energy based on user-specified EOS
  // currently, computing eps for classical ideal gas

  /* Computing v_j */
  const array<CCTK_REAL, 3> v_up = {pv.velx, pv.vely, pv.velz};
  const array<CCTK_REAL, 3> v_low =
      calc_vlow(v_up, g.gxx, g.gxy, g.gxz, g.gyy, g.gyz, g.gzz);

  const CCTK_REAL w_lorentz = calc_wlor(v_low, v_up);

  /* Computing beta_j */
  const array<CCTK_REAL, 3> beta_up = {shft.betax, shft.betay, shft.betaz};
  const array<CCTK_REAL, 3> beta_low =
      calc_vlow(beta_up, g.gxx, g.gxy, g.gxz, g.gyy, g.gyz, g.gzz);

  /* Computing B_j */
  const array<CCTK_REAL, 3> B_up = {pv.Bvecx, pv.Bvecy, pv.Bvecz};
  const array<CCTK_REAL, 3> B_low =
      calc_vlow(B_up, g.gxx, g.gxy, g.gxz, g.gyy, g.gyz, g.gzz);

  /* Computing b^t : this is b^0 * alp */
  const CCTK_REAL bst = w_lorentz * (B_low[0] * v_low[0] + B_low[1] * v_low[1] +
                                     B_low[2] * v_low[2]);

  /* Computing b^j */
  const CCTK_REAL bsx =
      (B_up[0] + bst * w_lorentz * (v_up[0] - beta_up[0] / lap.alp)) /
      w_lorentz;
  const CCTK_REAL bsy =
      (B_up[1] + bst * w_lorentz * (v_up[1] - beta_up[1] / lap.alp)) /
      w_lorentz;
  const CCTK_REAL bsz =
      (B_up[2] + bst * w_lorentz * (v_up[2] - beta_up[2] / lap.alp)) /
      w_lorentz;

  /* Computing b_j */
  const array<CCTK_REAL, 3> b_up = {bsx, bsy, bsz};
  array<CCTK_REAL, 3> b_low =
      calc_vlow(b_up, g.gxx, g.gxy, g.gxz, g.gyy, g.gyz, g.gzz);

  b_low[0] += beta_low[0] * bst / lap.alp;
  b_low[1] += beta_low[1] * bst / lap.alp;
  b_low[2] += beta_low[2] * bst / lap.alp;

  /* Computing b^mu b_mu */
  const CCTK_REAL bs2 = (B_up[0] * B_low[0] + B_up[1] * B_low[1] +
                         B_up[2] * B_low[2] + bst * bst) /
                        (w_lorentz * w_lorentz);

  // computing conserved from primitives
  cv.dens = sqrt_detg * pv.rho * w_lorentz;

  cv.momx =
      sqrt_detg * (w_lorentz * w_lorentz *
                       (pv.rho * (1 + pv.eps) + pv.press + bs2) * v_low[0] -
                   bst * b_low[0]);

  cv.momy =
      sqrt_detg * (w_lorentz * w_lorentz *
                       (pv.rho * (1 + pv.eps) + pv.press + bs2) * v_low[1] -
                   bst * b_low[1]);

  cv.momz =
      sqrt_detg * (w_lorentz * w_lorentz *
                       (pv.rho * (1 + pv.eps) + pv.press + bs2) * v_low[2] -
                   bst * b_low[2]);

  cv.tau = sqrt_detg * (w_lorentz * w_lorentz *
                            (pv.rho * (1 + pv.eps) + pv.press + bs2) -
                        (pv.press + 0.5 * bs2) - bst * bst) -
           cv.dens;

  cv.dBvecx = sqrt_detg * pv.Bvecx;
  cv.dBvecy = sqrt_detg * pv.Bvecy;
  cv.dBvecz = sqrt_detg * pv.Bvecz;
}

} // namespace AsterX
