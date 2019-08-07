#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop.hxx>

#include <cassert>
#include <cmath>
#include <iostream>

namespace HydroToyAMReX {
using namespace std;

constexpr int dim = 3;

////////////////////////////////////////////////////////////////////////////////

extern "C" void HydroToyAMReX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  CCTK_REAL x0[dim], dx[dim];
  for (int d = 0; d < dim; ++d) {
    x0[d] = CCTK_ORIGIN_SPACE(d);
    dx[d] = CCTK_DELTA_SPACE(d);
  }

  Loop::loop_all(cctkGH, [&](int i, int j, int k, int idx) {
    rho[idx] = 1.0;
    momx[idx] = 0.0;
    momy[idx] = 0.0;
    momz[idx] = 0.0;
    etot[idx] = 1.0;
  });
}

extern "C" void HydroToyAMReX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  constexpr int di = 1;
  const int dj = di * cctk_ash[0];
  const int dk = dj * cctk_ash[1];

  Loop::loop_int(cctkGH, [&](int i, int j, int k, int idx) {
    // dt rho + d_i (rho vel^i) = 0
    rho[idx] = rho_p[idx] - dt * ((momx_p[idx + di] - momx_p[idx]) / dx +
                                  (momy_p[idx + dj] - momy_p[idx]) / dy +
                                  (momz_p[idx + dk] - momz_p[idx]) / dz);

    //     // dt mom^i + d_j (mom^i vel^j) = 0
    //     momx[idx] = momx_p[idx] - dt * (
    // momx_p[idx] * momx_p[idx]
    //                                     );

    CCTK_REAL velx = momx[idx] / rho[idx];
    CCTK_REAL vely = momy[idx] / rho[idx];
    CCTK_REAL velz = momz[idx] / rho[idx];
  });
}

} // namespace HydroToyAMReX
