#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

extern "C" void HydroToyCarpetX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Evolve;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  // Transport
  // dt rho + d_i (rho vel^i) = 0
  // dt mom_j + d_i (mom_j vel^i) = 0
  // dt etot + d_i (etot vel^i) = 0

  const CCTK_REAL dt_dx = dt / dx;
  const CCTK_REAL dt_dy = dt / dy;
  const CCTK_REAL dt_dz = dt / dz;

#if 0
    // CPU

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      auto calcupdate{[&](auto &fx, auto &fy, auto &fz) {
        return dt_dx * (fx(p.I + p.DI(0)) - fx(p.I)) +
               dt_dy * (fy(p.I + p.DI(1)) - fy(p.I)) +
               dt_dz * (fz(p.I + p.DI(2)) - fz(p.I));
      }};

      rho_(p.I) = rho_p_(p.I) - calcupdate(fxrho_, fyrho_, fzrho_);

      // momx_(p.I) = momx_p_(p.I) - calcupdate(fxmomx_, fymomx_, fzmomx_);
      // momy_(p.I) = momy_p_(p.I) - calcupdate(fxmomy_, fymomy_, fzmomy_);
      // momz_(p.I) = momz_p_(p.I) - calcupdate(fxmomz_, fymomz_, fzmomz_);

      // etot_(p.I) = etot_p_(p.I) - calcupdate(fxetot_, fyetot_, fzetot_);
    });
#else
  // CPU or GPU
  // TODO: &p?
  //  const auto calcupdate =
  //      [=] CCTK_DEVICE CCTK_HOST (CCTK_REAL fx, CCTK_REAL fy, CCTK_REAL fz, &p)
  //      CCTK_ATTRIBUTE_ALWAYS_INLINE {
  //	  return dt_dx * (fx(p.idx + p.di) - fx(p.idx)) +
  //             dt_dy * (fy(p.idx + p.dj) - fy(p.idx)) +
  //             dt_dz * (fz(p.idx + p.dk) - fz(p.idx));
  //      };

  printf("Before loop\n");

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);
  griddesc.loop_int_device<1, 1, 1>(
      nghostzones, [=] CCTK_DEVICE  CCTK_HOST (const Loop::PointDesc &p)
                       CCTK_ATTRIBUTE_ALWAYS_INLINE{
                         const auto px = griddesc.point_desc<0, 1, 1>(p);
                         const auto py = griddesc.point_desc<1, 0, 1>(p);
                         const auto pz = griddesc.point_desc<1, 1, 0>(p);
                         CCTK_REAL calcupdate =
                             dt_dx * (fxrho[px.idx + px.di] - fxrho[px.idx]) +
                             dt_dy * (fyrho[py.idx + py.dj] - fyrho[py.idx]) +
                             dt_dz * (fzrho[pz.idx + pz.dk] - fzrho[pz.idx]);

                         rho[p.idx] = rho_p[p.idx] - calcupdate;
                         //        rho[p.idx] = rho_p[p.idx] - calcupdate(fxrho,
                         //        fyrho, fzrho, &p);

                         //        psi[p.idx] =
                         //            psi_p[p.idx] +
                         //            dt * (ddx_phi + ddy_phi + ddz_phi -
                         //            pow(mass, 2) * phi_p[p.idx] +
                         //                  4 * M_PI * central_potential(t,
                         //                  p.x, p.y, p.z));
                         //        phi[p.idx] = phi_p[p.idx] + dt * psi[p.idx];
                       });
#endif
}

} // namespace HydroToyCarpetX
