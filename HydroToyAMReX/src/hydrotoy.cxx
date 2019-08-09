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

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    rho[p.idx] = 1.0;
    momx[p.idx] = 0.0;
    momy[p.idx] = 0.0;
    momz[p.idx] = 0.0;
    etot[p.idx] = 1.0;
  });
}

extern "C" void HydroToyAMReX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;

  constexpr int di = 1;
  const int dj = di * cctk_ash[0];
  const int dk = dj * cctk_ash[1];

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    // dt rho + d_i (rho vel^i) = 0
    rho[p.idx] =
        rho_p[p.idx] - dt * ((momx_p[p.idx + di] - momx_p[p.idx]) / p.dx +
                             (momy_p[p.idx + dj] - momy_p[p.idx]) / p.dy +
                             (momz_p[p.idx + dk] - momz_p[p.idx]) / p.dz);

    //     // dt mom^i + d_j (mom^i vel^j) = 0
    //     momx[p.idx] = momx_p[p.idx] - dt * (
    // momx_p[p.idx] * momx_p[p.idx]
    //                                     );

    CCTK_REAL velx = momx[p.idx] / rho[p.idx];
    CCTK_REAL vely = momy[p.idx] / rho[p.idx];
    CCTK_REAL velz = momz[p.idx] / rho[p.idx];
  });
}

} // namespace HydroToyAMReX
