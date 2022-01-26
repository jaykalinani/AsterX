#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

struct metric {
  CCTK_REAL gxx, gxy, gxz, gyy, gyz, gzz;
};

struct prim {
  CCTK_REAL rho;
  CCTK_REAL velx, vely, velz;
  CCTK_REAL eps, press;
};

struct cons {
  CCTK_REAL dens;
  CCTK_REAL momx, momy, momz;
  CCTK_REAL tau;
};

CCTK_DEVICE CCTK_HOST void prim2con(const metric &g, const prim &pv, cons &cv) {

  // determinant of spatial metric
  const CCTK_REAL detg = -g.gxz * g.gxz * g.gyy + 2.0 * g.gxy * g.gxz * g.gyz -
                         g.gxx * g.gyz * g.gyz - g.gxy * g.gxy * g.gzz +
                         g.gxx * g.gyy * g.gzz;
  const CCTK_REAL sqrt_detg = sqrt(detg);

  // TODO: compute specific internal energy based on user-specified EOS
  // currently, computing eps for classical ideal gas

  // v_j
  const CCTK_REAL vlowx = g.gxx * pv.velx + g.gxy * pv.vely + g.gxz * pv.velz;
  const CCTK_REAL vlowy = g.gxy * pv.velx + g.gyy * pv.vely + g.gyz * pv.velz;
  const CCTK_REAL vlowz = g.gxz * pv.velx + g.gyz * pv.vely + g.gzz * pv.velz;

  // w_lorentz
  const CCTK_REAL w_lorentz =
      1.0 /
      sqrt(1 - (g.gxx * pv.velx * pv.velx + g.gyy * pv.vely * pv.vely +
                g.gzz * pv.velz * pv.velz + 2 * g.gxy * pv.velx * pv.vely +
                2 * g.gxz * pv.velx * pv.velz + 2 * g.gyz * pv.vely * pv.velz));

  // computing conserved from primitives
  cv.dens = sqrt_detg * pv.rho * w_lorentz;

  cv.momx =
      sqrt_detg * pv.rho * w_lorentz * (1 + pv.eps + pv.press / pv.rho) * vlowx;

  cv.momy =
      sqrt_detg * pv.rho * w_lorentz * (1 + pv.eps + pv.press / pv.rho) * vlowy;

  cv.momz =
      sqrt_detg * pv.rho * w_lorentz * (1 + pv.eps + pv.press / pv.rho) * vlowz;

  cv.tau = sqrt_detg * pv.rho * w_lorentz *
               ((1 + pv.eps + pv.press / pv.rho) * w_lorentz - 1) -
           pv.press;
}

} // namespace GRHydroToyGPU
