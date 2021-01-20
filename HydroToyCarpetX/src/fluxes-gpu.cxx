//#include <loop.hxx>

//#include <vectors.h>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>
#include <loop_device.hxx>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>

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

  // frho^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // fetot^i = (etot + press) vel^i

  const auto calcflux =
          [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
              CCTK_REAL var_p, CCTK_REAL var_m, CCTK_REAL flux_p, CCTK_REAL flux_m) {

            CCTK_REAL lambda_m = 1.0;
            CCTK_REAL lambda_p = -1.0;
//            CCTK_REAL var_m = u[p.idx-p.di];
//            CCTK_REAL var_p = u[p.idx];
//            CCTK_REAL flux_m = f[p1.idx-p1.di];
//            CCTK_REAL flux_p = f[p1.idx];
            CCTK_REAL llf = 0.5 *((flux_m + flux_p) - fmax(fabs(lambda_m), fabs(lambda_p)) * (var_p - var_m));
//            printf("llf = %g, dAx = %g\n",llf,dAx);
            return dAx * llf;
          };

//  printf("Before loop\n");

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};

  const Loop::GridDescBaseDevice griddesc(cctkGH);
  
  griddesc.loop_all_device<1, 1, 1>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {
      const auto p1 = griddesc.point_desc<0, 1, 1>(p);
//      const auto p2 = griddesc.point_desc<1, 1, 1>(p);
//      printf("p: x = %g, y = %g, z = %g\n", p.x, p.y, p.z);
//      printf("p1: x = %g, y = %g, z = %g\n", p1.x, p1.y, p1.z);
//      printf("p: i = %d, j = %d, k = %d\n", p.i, p.j, p.k);
//      printf("p.idx = %d, p1.idx = %d\n", p.idx, p1.idx);
//      printf("p1: i = %d, j = %d, k = %d\n", p1.i, p1.j, p1.k);
//      assert(p1.idx == p.idx);
//      fxrho[p1.idx] = calcflux(rho, rho[p.idx-p.di], p, p1);
//      fxrho[p1.idx] = calcflux(rho[p.idx], rho[p.idx-p.di], p, p1);
      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("p: x = %g, y = %g, z = %g\n", p.x, p.y, p.z);
        printf("p1: x = %g, y = %g, z = %g\n", p1.x, p1.y, p1.z);
//        printf("p2: x = %g, y = %g, z = %g\n", p2.x, p2.y, p2.z);
        printf("p: i = %d, j = %d, k = %d\n", p.i, p.j, p.k);
        printf("p1: i = %d, j = %d, k = %d\n", p1.i, p1.j, p1.k);
//        printf("p2: i = %d, j = %d, k = %d\n", p2.i, p2.j, p2.k);
        printf("rho_p = %g, rho_m = %g\n", rho[p.idx], rho[p.idx-p.di]);
        printf("momx_p = %g, momx_m = %g\n", momx[p.idx], momx[p.idx-p.di]);
        printf("momy_p = %g, momy_m = %g\n", momy[p.idx], momy[p.idx-p.di]);
        printf("momz_p = %g, momz_m = %g\n", momz[p.idx], momz[p.idx-p.di]);
        printf("etot_p = %g, etot_m = %g\n", etot[p.idx], etot[p.idx-p.di]);
        printf("fxrho_p = %f, fxrho_m =%f\n",(rho[p.idx] * velx[p.idx]),(rho[p.idx-p.di] * velx[p.idx-p.di]));
        printf("fxmomx_p = %f, fxmomx_m =%f\n",(momx[p.idx] * velx[p.idx] + press[p.idx],momx[p.idx-p.di] * velx[p.idx-p.di] + press[p.idx-p.di]));
        printf("fxmomy_p = %f, fxmomy_m =%f\n",(momy[p.idx] * velx[p.idx] + press[p.idx],momx[p.idx-p.di] * velx[p.idx-p.di] + press[p.idx-p.di]));
        printf("fxmomz_p = %f, fxmomz_m =%f\n",(momz[p.idx] * velx[p.idx] + press[p.idx],momx[p.idx-p.di] * velx[p.idx-p.di] + press[p.idx-p.di]));
        printf("fxetot_p = %f, fxetot_m =%f\n",(etot[p.idx] + press[p.idx]) * velx[p.idx], (etot[p.idx-p.di] + press[p.idx-p.di]) * velx[p.idx-p.di]);
      }

      fxrho[p1.idx] = calcflux(rho[p.idx], rho[p.idx-p.di], rho[p.idx] * velx[p.idx], rho[p.idx-p.di] * velx[p.idx-p.di]);
      fxmomx[p1.idx] = calcflux(momx[p.idx], momx[p.idx-p.di],momx[p.idx] * velx[p.idx] + press[p.idx],momx[p.idx-p.di] * velx[p.idx-p.di] + press[p.idx-p.di]);
      fxmomy[p1.idx] = calcflux(momy[p.idx], momy[p.idx-p.di],momy[p.idx] * velx[p.idx] + press[p.idx],momy[p.idx-p.di] * velx[p.idx-p.di] + press[p.idx-p.di]);
      fxmomz[p1.idx] = calcflux(momz[p.idx], momz[p.idx-p.di],momz[p.idx] * velx[p.idx] + press[p.idx],momz[p.idx-p.di] * velx[p.idx-p.di] + press[p.idx-p.di]);
      fxetot[p1.idx] = calcflux(etot[p.idx], etot[p.idx-p.di], (etot[p.idx] + press[p.idx]) * velx[p.idx], (etot[p.idx-p.di] + press[p.idx-p.di]) * velx[p.idx-p.di]);

      fxrho[p.idx] = 0.0;
      fxmomx[p.idx] = 0.0;
      fxmomy[p.idx] = 0.0;
      fxmomz[p.idx] = 0.0;
      fxetot[p.idx] = 0.0;
 
      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("fxrho = %f\n",fxrho[p1.idx]);
        printf("fxmomx = %f\n",fxmomx[p1.idx]);
        printf("fxmomy = %f\n",fxmomy[p1.idx]);
        printf("fxmomz = %f\n",fxmomz[p1.idx]);
        printf("fxetot = %f\n",fxetot[p1.idx]);
      }

  });

  griddesc.loop_all_device<1, 1, 1>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {
      const auto p1 = griddesc.point_desc<1, 0, 1>(p);
      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("p: x = %g, y = %g, z = %g\n", p.x, p.y, p.z);
        printf("p1: x = %g, y = %g, z = %g\n", p1.x, p1.y, p1.z);
        printf("p: i = %d, j = %d, k = %d\n", p.i, p.j, p.k);
        printf("p1: i = %d, j = %d, k = %d\n", p1.i, p1.j, p1.k);
        printf("rho_p = %g, rho_m = %g\n", rho[p.idx], rho[p.idx-p.dj]);
        printf("momx_p = %g, momx_m = %g\n", momx[p.idx], momx[p.idx-p.di]);
        printf("momy_p = %g, momy_m = %g\n", momy[p.idx], momy[p.idx-p.di]);
        printf("momz_p = %g, momz_m = %g\n", momz[p.idx], momz[p.idx-p.di]);
        printf("etot_p = %g, etot_m = %g\n", etot[p.idx], etot[p.idx-p.di]);
        printf("fyrho_p = %f, fyrho_m =%f\n",(rho[p.idx] * vely[p.idx]),(rho[p.idx-p.dj] * vely[p.idx-p.dj]));
        printf("fymomx_p = %f, fymomx_m =%f\n",(momx[p.idx] * vely[p.idx] + press[p.idx],momx[p.idx-p.dj] * vely[p.idx-p.dj] + press[p.idx-p.dj]));
        printf("fymomy_p = %f, fymomy_m =%f\n",(momy[p.idx] * vely[p.idx] + press[p.idx],momx[p.idx-p.dj] * vely[p.idx-p.dj] + press[p.idx-p.dj]));
        printf("fymomz_p = %f, fymomz_m =%f\n",(momz[p.idx] * vely[p.idx] + press[p.idx],momx[p.idx-p.dj] * vely[p.idx-p.dj] + press[p.idx-p.dj]));
        printf("fyetot_p = %f, fyetot_m =%f\n",(etot[p.idx] + press[p.idx]) * vely[p.idx], (etot[p.idx-p.dj] + press[p.idx-p.dj]) * vely[p.idx-p.dj]);
      }

      fyrho[p1.idx] = calcflux(rho[p.idx], rho[p.idx-p.dj], rho[p.idx] * vely[p.idx], rho[p.idx-p.dj] * vely[p.idx-p.dj]);
      fymomx[p1.idx] = calcflux(momx[p.idx], momx[p.idx-p.dj],momx[p.idx] * vely[p.idx] + press[p.idx],momx[p.idx-p.dj] * vely[p.idx-p.dj] + press[p.idx-p.dj]);
      fymomy[p1.idx] = calcflux(momy[p.idx], momy[p.idx-p.dj],momy[p.idx] * vely[p.idx] + press[p.idx],momy[p.idx-p.dj] * vely[p.idx-p.dj] + press[p.idx-p.dj]);
      fymomz[p1.idx] = calcflux(momz[p.idx], momz[p.idx-p.dj],momz[p.idx] * vely[p.idx] + press[p.idx],momz[p.idx-p.dj] * vely[p.idx-p.dj] + press[p.idx-p.dj]);
      fyetot[p1.idx] = calcflux(etot[p.idx], etot[p.idx-p.dj], (etot[p.idx] + press[p.idx]) * vely[p.idx], (etot[p.idx-p.dj] + press[p.idx-p.dj]) * vely[p.idx-p.dj]);

      fyrho[p.idx] = 0.0;
      fymomx[p.idx] = 0.0;
      fymomy[p.idx] = 0.0;
      fymomz[p.idx] = 0.0;
      fyetot[p.idx] = 0.0;
 
      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("fyrho = %f\n",fyrho[p1.idx]);
        printf("fymomx = %f\n",fymomx[p1.idx]);
        printf("fymomy = %f\n",fymomy[p1.idx]);
        printf("fymomz = %f\n",fymomz[p1.idx]);
        printf("fyetot = %f\n",fyetot[p1.idx]);
      }

  });

  griddesc.loop_all_device<1, 1, 1>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {
      const auto p1 = griddesc.point_desc<1, 1, 0>(p);
      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("p: x = %g, y = %g, z = %g\n", p.x, p.y, p.z);
        printf("p1: x = %g, y = %g, z = %g\n", p1.x, p1.y, p1.z);
        printf("p: i = %d, j = %d, k = %d\n", p.i, p.j, p.k);
        printf("p1: i = %d, j = %d, k = %d\n", p1.i, p1.j, p1.k);
        printf("rho_p = %g, rho_m = %g\n", rho[p.idx], rho[p.idx-p.dk]);
        printf("momx_p = %g, momx_m = %g\n", momx[p.idx], momx[p.idx-p.dk]);
        printf("momy_p = %g, momy_m = %g\n", momy[p.idx], momy[p.idx-p.dk]);
        printf("momz_p = %g, momz_m = %g\n", momz[p.idx], momz[p.idx-p.dk]);
        printf("etot_p = %g, etot_m = %g\n", etot[p.idx], etot[p.idx-p.dk]);
        printf("fzrho_p = %f, fzrho_m =%f\n",(rho[p.idx] * velz[p.idx]),(rho[p.idx-p.dk] * velz[p.idx-p.dk]));
        printf("fzmomx_p = %f, fzmomx_m =%f\n",(momx[p.idx] * velz[p.idx] + press[p.idx],momx[p.idx-p.dk] * velz[p.idx-p.dk] + press[p.idx-p.dk]));
        printf("fzmomy_p = %f, fzmomy_m =%f\n",(momy[p.idx] * velz[p.idx] + press[p.idx],momx[p.idx-p.dk] * velz[p.idx-p.dk] + press[p.idx-p.dk]));
        printf("fzmomz_p = %f, fzmomz_m =%f\n",(momz[p.idx] * velz[p.idx] + press[p.idx],momx[p.idx-p.dk] * velz[p.idx-p.dk] + press[p.idx-p.dk]));
        printf("fzetot_p = %f, fzetot_m =%f\n",(etot[p.idx] + press[p.idx]) * velz[p.idx], (etot[p.idx-p.dk] + press[p.idx-p.dk]) * velz[p.idx-p.dk]);
      }

      fzrho[p1.idx] = calcflux(rho[p.idx], rho[p.idx-p.dk], rho[p.idx] * velz[p.idx], rho[p.idx-p.dk] * velz[p.idx-p.dk]);
      fzmomx[p1.idx] = calcflux(momx[p.idx], momx[p.idx-p.dk],momx[p.idx] * velz[p.idx] + press[p.idx],momx[p.idx-p.dk] * velz[p.idx-p.dk] + press[p.idx-p.dk]);
      fzmomy[p1.idx] = calcflux(momy[p.idx], momy[p.idx-p.dk],momy[p.idx] * velz[p.idx] + press[p.idx],momy[p.idx-p.dk] * velz[p.idx-p.dk] + press[p.idx-p.dk]);
      fzmomz[p1.idx] = calcflux(momz[p.idx], momz[p.idx-p.dk],momz[p.idx] * velz[p.idx] + press[p.idx],momz[p.idx-p.dk] * velz[p.idx-p.dk] + press[p.idx-p.dk]);
      fzetot[p1.idx] = calcflux(etot[p.idx], etot[p.idx-p.dk], (etot[p.idx] + press[p.idx]) * velz[p.idx], (etot[p.idx-p.dk] + press[p.idx-p.dk]) * velz[p.idx-p.dk]);

      fzrho[p.idx] = 0.0;
      fzmomx[p.idx] = 0.0;
      fzmomy[p.idx] = 0.0;
      fzmomz[p.idx] = 0.0;
      fzetot[p.idx] = 0.0;
 
      if (p1.i == 3 && p1.j == 4 && p1.k == 5) {
        printf("fzrho = %f\n",fzrho[p1.idx]);
        printf("fzmomx = %f\n",fzmomx[p1.idx]);
        printf("fzmomy = %f\n",fzmomy[p1.idx]);
        printf("fzmomz = %f\n",fzmomz[p1.idx]);
        printf("fzetot = %f\n",fzetot[p1.idx]);
      }

  });

//  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
//    auto calcflux=[&](auto &u, auto f) {
//      auto I_m = p.I - p.DI(1);
//      auto I_p = p.I;
//      auto lambda_m = 1.0;
//      auto lambda_p = -1.0;
//      auto var_m = u(I_m);
//      auto var_p = u(I_p);
//      auto flux_m = f(I_m);
//      auto flux_p = f(I_p);
//      return dAy * llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
//    };
//
//    fyrho_(p.I) = calcflux(rho_, [&](auto I) { return rho_(I) * vely_(I); });
//
//    fymomx_(p.I) = calcflux(momx_, [&](auto I) { return momx_(I) * vely_(I); });
//    fymomy_(p.I) = calcflux(
//        momy_, [&](auto I) { return momy_(I) * vely_(I) + press_(I); });
//    fymomz_(p.I) = calcflux(momz_, [&](auto I) { return momz_(I) * vely_(I); });
//
//    fyetot_(p.I) = calcflux(
//        etot_, [&](auto I) { return (etot_(I) + press_(I)) * vely_(I); });
//  });
//
//  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
//    auto calcflux=[&](auto &u, auto f) {
//      auto I_m = p.I - p.DI(2);
//      auto I_p = p.I;
//      auto lambda_m = 1.0;
//      auto lambda_p = -1.0;
//      auto var_m = u(I_m);
//      auto var_p = u(I_p);
//      auto flux_m = f(I_m);
//      auto flux_p = f(I_p);
//      return dAz * llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
//    };
//
//    fzrho_(p.I) = calcflux(rho_, [&](auto I) { return rho_(I) * velz_(I); });
//
//    fzmomx_(p.I) = calcflux(momx_, [&](auto I) { return momx_(I) * velz_(I); });
//    fzmomy_(p.I) = calcflux(momy_, [&](auto I) { return momy_(I) * velz_(I); });
//    fzmomz_(p.I) = calcflux(
//        momz_, [&](auto I) { return momz_(I) * velz_(I) + press_(I); });
//
//    fzetot_(p.I) = calcflux(
//        etot_, [&](auto I) { return (etot_(I) + press_(I)) * velz_(I); });
//  });
 }
} // namespace HydroToyCarpetX
