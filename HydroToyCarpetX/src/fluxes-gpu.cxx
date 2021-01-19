#include <loop.hxx>

#include <vectors.h>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

extern "C" void HydroToyCarpetX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dAx = dt * dy * dz;
  const CCTK_REAL dAy = dt * dx * dz;
  const CCTK_REAL dAz = dt * dx * dy;

  // Fluxes
  // frho^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // fetot^i = (etot + press) vel^i

  // For refluxing, AMReX expects fluxes to be integrated over the
  // faces both in space and time.
  const auto calcflux =
      [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
          CCTK_REAL var_p, CCTK_REAL var_m, CCTK_REAL f_p, CCTK_REAL f_m) {
            CCTK_REAL lambda_m = 1.0;
            CCTK_REAL lambda_p  -1.0;
            llf = 0.5 * ((flux_m + flux_p) - fmax(fabs(lambda_m), fabs(lambda_p)) * (var_p - var_m));
      	    return dAx * llf;
      };

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);
  griddesc.loop_all_device<1, 1, 1>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {
        const auto p1 = griddesc.point_desc<0, 1, 1>(p);
        fxrho[p1.idx] = calcflux(rho[p.idx],rho[p.idx-p.di],rho[p.idx]*velx[p.idx],rho[p.idx-p.di]*velx[p.idx-p.di]);
//        fxrho[p.idx] = calcflux(rho[p.idx],rho[p.idx-p.di],rho[p.idx]*velx[p.idx],rho[p.idx-p.di]*velx[p.idx-p.di]); if I set this I don't get NaNs.
    };
  }
} // namespace HydroToyCarpetX
