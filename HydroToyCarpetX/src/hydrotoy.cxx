#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <iostream>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

////////////////////////////////////////////////////////////////////////////////

template <typename T> inline T fmax5(T x0, T x1, T x2, T x3, T x4) {
  T r0 = fmax(x0, x1);
  T r1 = fmax(x2, x3);
  T r2 = fmax(x4, r0);
  T r3 = fmax(r1, r2);
  return r3;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void HydroToyCarpetX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D1<CCTK_REAL> rho_(cctkGH, {1, 1, 1}, {1, 1, 1}, rho);
  const Loop::GF3D1<CCTK_REAL> momx_(cctkGH, {1, 1, 1}, {1, 1, 1}, momx);
  const Loop::GF3D1<CCTK_REAL> momy_(cctkGH, {1, 1, 1}, {1, 1, 1}, momy);
  const Loop::GF3D1<CCTK_REAL> momz_(cctkGH, {1, 1, 1}, {1, 1, 1}, momz);
  const Loop::GF3D1<CCTK_REAL> etot_(cctkGH, {1, 1, 1}, {1, 1, 1}, etot);

  if (CCTK_EQUALS(setup, "equilibrium")) {

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      rho_(p.I) = 1.0;
      momx_(p.I) = 0.0;
      momy_(p.I) = 0.0;
      momz_(p.I) = 0.0;
      etot_(p.I) = 1.0;
    });

  } else if (CCTK_EQUALS(setup, "sound wave")) {

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      rho_(p.I) = 1.0;
      momx_(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
      momy_(p.I) = 0.0;
      momz_(p.I) = 0.0;
      etot_(p.I) = 1.0; // should add kinetic energy here
    });

  } else if (CCTK_EQUALS(setup, "shock tube")) {

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      if (p.x <= 0.0) {
        rho_(p.I) = 2.0;
        momx_(p.I) = 0.0;
        momy_(p.I) = 0.0;
        momz_(p.I) = 0.0;
        etot_(p.I) = 2.0;
      } else {
        rho_(p.I) = 1.0;
        momx_(p.I) = 0.0;
        momy_(p.I) = 0.0;
        momz_(p.I) = 0.0;
        etot_(p.I) = 1.0;
      }
    });

  } else if (CCTK_EQUALS(setup, "spherical shock")) {

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      CCTK_REAL r2 = pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2);
      if (r2 <= pow(shock_radius, 2)) {
        rho_(p.I) = 2.0;
        momx_(p.I) = 0.0;
        momy_(p.I) = 0.0;
        momz_(p.I) = 0.0;
        etot_(p.I) = 2.0;
      } else {
        rho_(p.I) = 1.0;
        momx_(p.I) = 0.0;
        momy_(p.I) = 0.0;
        momz_(p.I) = 0.0;
        etot_(p.I) = 1.0;
      }
    });

  } else {
    assert(0);
  }
}

extern "C" void HydroToyCarpetX_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Boundaries;

  // do nothing

  Loop::loop_bnd<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    assert(false); // This should not be executed
  });
}

extern "C" void HydroToyCarpetX_CopyConserved(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_CopyConserved;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D1<const CCTK_REAL> rho_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                            rho_p);
  const Loop::GF3D1<const CCTK_REAL> momx_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momx_p);
  const Loop::GF3D1<const CCTK_REAL> momy_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momy_p);
  const Loop::GF3D1<const CCTK_REAL> momz_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momz_p);
  const Loop::GF3D1<const CCTK_REAL> etot_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             etot_p);

  const Loop::GF3D1<CCTK_REAL> rho_(cctkGH, {1, 1, 1}, {1, 1, 1}, rho);
  const Loop::GF3D1<CCTK_REAL> momx_(cctkGH, {1, 1, 1}, {1, 1, 1}, momx);
  const Loop::GF3D1<CCTK_REAL> momy_(cctkGH, {1, 1, 1}, {1, 1, 1}, momy);
  const Loop::GF3D1<CCTK_REAL> momz_(cctkGH, {1, 1, 1}, {1, 1, 1}, momz);
  const Loop::GF3D1<CCTK_REAL> etot_(cctkGH, {1, 1, 1}, {1, 1, 1}, etot);

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    rho_(p.I) = rho_p_(p.I);
    momx_(p.I) = momx_p_(p.I);
    momy_(p.I) = momy_p_(p.I);
    momz_(p.I) = momz_p_(p.I);
    etot_(p.I) = etot_p_(p.I);
  });
}

#if 0
extern "C" void HydroToyCarpetX_Pressure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Pressure;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D1<const CCTK_REAL> rho_(cctkGH, {1, 1, 1}, {1, 1, 1}, rho);
  const Loop::GF3D1<const CCTK_REAL> momx_(cctkGH, {1, 1, 1}, {1, 1, 1}, momx);
  const Loop::GF3D1<const CCTK_REAL> momy_(cctkGH, {1, 1, 1}, {1, 1, 1}, momy);
  const Loop::GF3D1<const CCTK_REAL> momz_(cctkGH, {1, 1, 1}, {1, 1, 1}, momz);
  const Loop::GF3D1<const CCTK_REAL> etot_(cctkGH, {1, 1, 1}, {1, 1, 1}, etot);

  const Loop::GF3D1<CCTK_REAL> press_(cctkGH, {1, 1, 1}, {1, 1, 1}, press);
  const Loop::GF3D1<CCTK_REAL> velx_(cctkGH, {1, 1, 1}, {1, 1, 1}, velx);
  const Loop::GF3D1<CCTK_REAL> vely_(cctkGH, {1, 1, 1}, {1, 1, 1}, vely);
  const Loop::GF3D1<CCTK_REAL> velz_(cctkGH, {1, 1, 1}, {1, 1, 1}, velz);
  const Loop::GF3D1<CCTK_REAL> eint_(cctkGH, {1, 1, 1}, {1, 1, 1}, eint);

  // Equation of state: p = (gamma - 1) e

  // vel^j = delta^j_i mom_i / rho

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    CCTK_REAL ekin =
        sqrt(pow(momx_(p.I), 2) + pow(momy_(p.I), 2) + pow(momz_(p.I), 2)) /
        (2 * rho_(p.I));
    eint_(p.I) = etot_(p.I) - ekin;

    press_(p.I) = (gamma - 1) * eint_(p.I);

    velx_(p.I) = momx_(p.I) / rho_(p.I);
    vely_(p.I) = momy_(p.I) / rho_(p.I);
    velz_(p.I) = momz_(p.I) / rho_(p.I);
  });
}
#endif

// Minmod slope limiter
template <typename T> T minmod(T a, T b) {
  // Same sign: return value with smaller magnitude
  // Different signs: return zero
  return (copysign(T(0.5), a) + copysign(T(0.5), b)) * fmin(fabs(a), fabs(b));
}

// Van Leer slope limiter
template <typename T> T vanleer(T a, T b) {
  return ifthen(a * b > 0, 2 * a * b / (a + b), T(0));
}

// MC2 slope limiter
template <typename T> T mc2(T a, T b) {
  return (copysign(T(0.25), a) + copysign(T(0.25), b)) *
         fmin(4 * fmin(fabs(a), fabs(b)), fabs(a) + fabs(b));
}

// Superbee slope limiter
template <typename T> T superbee(T a, T b) {
  return (copysign(T(0.5), a) + copysign(T(0.5), b)) *
         fmax(fmin(2 * fabs(a), fabs(b)), fmin(fabs(a), 2 * fabs(b)));
}

// LLF (local Lax-Friedrichs) Riemann solver
template <typename T>
T llf(T lambda_m, T lambda_p, T var_m, T var_p, T flux_m, T flux_p) {
  return T(0.5) * ((flux_m + flux_p) -
                   fmax(fabs(lambda_m), fabs(lambda_p)) * (var_p - var_m));
}

// HLLE (Harten, Lax, van Leer, Einfeldt) Riemann solver
template <typename T>
T hlle(T lambda_m, T lambda_p, T var_m, T var_p, T flux_m, T flux_p) {
  // var_m and var_p are probably swapped
  assert(0);
  return (lambda_p * flux_p - lambda_m * flux_m +
          lambda_p * lambda_m * (var_m - var_p)) /
         (lambda_p - lambda_m);
}

#if 0
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

  const Loop::GF3D1<const CCTK_REAL> rho_(cctkGH, {1, 1, 1}, {1, 1, 1}, rho);
  const Loop::GF3D1<const CCTK_REAL> momx_(cctkGH, {1, 1, 1}, {1, 1, 1}, momx);
  const Loop::GF3D1<const CCTK_REAL> momy_(cctkGH, {1, 1, 1}, {1, 1, 1}, momy);
  const Loop::GF3D1<const CCTK_REAL> momz_(cctkGH, {1, 1, 1}, {1, 1, 1}, momz);
  const Loop::GF3D1<const CCTK_REAL> etot_(cctkGH, {1, 1, 1}, {1, 1, 1}, etot);

  const Loop::GF3D1<const CCTK_REAL> press_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                            press);
  const Loop::GF3D1<const CCTK_REAL> velx_(cctkGH, {1, 1, 1}, {1, 1, 1}, velx);
  const Loop::GF3D1<const CCTK_REAL> vely_(cctkGH, {1, 1, 1}, {1, 1, 1}, vely);
  const Loop::GF3D1<const CCTK_REAL> velz_(cctkGH, {1, 1, 1}, {1, 1, 1}, velz);

  const Loop::GF3D1<CCTK_REAL> fxrho_(cctkGH, {0, 1, 1}, {0, 0, 0}, fxrho);
  const Loop::GF3D1<CCTK_REAL> fxmomx_(cctkGH, {0, 1, 1}, {0, 0, 0}, fxmomx);
  const Loop::GF3D1<CCTK_REAL> fxmomy_(cctkGH, {0, 1, 1}, {0, 0, 0}, fxmomy);
  const Loop::GF3D1<CCTK_REAL> fxmomz_(cctkGH, {0, 1, 1}, {0, 0, 0}, fxmomz);
  const Loop::GF3D1<CCTK_REAL> fxetot_(cctkGH, {0, 1, 1}, {0, 0, 0}, fxetot);

  const Loop::GF3D1<CCTK_REAL> fyrho_(cctkGH, {1, 0, 1}, {0, 0, 0}, fyrho);
  const Loop::GF3D1<CCTK_REAL> fymomx_(cctkGH, {1, 0, 1}, {0, 0, 0}, fymomx);
  const Loop::GF3D1<CCTK_REAL> fymomy_(cctkGH, {1, 0, 1}, {0, 0, 0}, fymomy);
  const Loop::GF3D1<CCTK_REAL> fymomz_(cctkGH, {1, 0, 1}, {0, 0, 0}, fymomz);
  const Loop::GF3D1<CCTK_REAL> fyetot_(cctkGH, {1, 0, 1}, {0, 0, 0}, fyetot);

  const Loop::GF3D1<CCTK_REAL> fzrho_(cctkGH, {1, 1, 0}, {0, 0, 0}, fzrho);
  const Loop::GF3D1<CCTK_REAL> fzmomx_(cctkGH, {1, 1, 0}, {0, 0, 0}, fzmomx);
  const Loop::GF3D1<CCTK_REAL> fzmomy_(cctkGH, {1, 1, 0}, {0, 0, 0}, fzmomy);
  const Loop::GF3D1<CCTK_REAL> fzmomz_(cctkGH, {1, 1, 0}, {0, 0, 0}, fzmomz);
  const Loop::GF3D1<CCTK_REAL> fzetot_(cctkGH, {1, 1, 0}, {0, 0, 0}, fzetot);

  // frho^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // fetot^i = (etot + press) vel^i

  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcflux=[&](auto &u, auto f) {
      auto I_m = p.I - p.DI(0);
      auto I_p = p.I;
      auto lambda_m = 1.0;
      auto lambda_p = -1.0;
      auto var_m = u(I_m);
      auto var_p = u(I_p);
      auto flux_m = f(I_m);
      auto flux_p = f(I_p);
      return dAx * llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
    };

    fxrho_(p.I) = calcflux(rho_, [&](auto I) { return rho_(I) * velx_(I); });

    fxmomx_(p.I) = calcflux(
        momx_, [&](auto I) { return momx_(I) * velx_(I) + press_(I); });
    fxmomy_(p.I) = calcflux(momy_, [&](auto I) { return momy_(I) * velx_(I); });
    fxmomz_(p.I) = calcflux(momz_, [&](auto I) { return momz_(I) * velx_(I); });

    fxetot_(p.I) = calcflux(
        etot_, [&](auto I) { return (etot_(I) + press_(I)) * velx_(I); });
  });

  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcflux=[&](auto &u, auto f) {
      auto I_m = p.I - p.DI(1);
      auto I_p = p.I;
      auto lambda_m = 1.0;
      auto lambda_p = -1.0;
      auto var_m = u(I_m);
      auto var_p = u(I_p);
      auto flux_m = f(I_m);
      auto flux_p = f(I_p);
      return dAy * llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
    };

    fyrho_(p.I) = calcflux(rho_, [&](auto I) { return rho_(I) * vely_(I); });

    fymomx_(p.I) = calcflux(momx_, [&](auto I) { return momx_(I) * vely_(I); });
    fymomy_(p.I) = calcflux(
        momy_, [&](auto I) { return momy_(I) * vely_(I) + press_(I); });
    fymomz_(p.I) = calcflux(momz_, [&](auto I) { return momz_(I) * vely_(I); });

    fyetot_(p.I) = calcflux(
        etot_, [&](auto I) { return (etot_(I) + press_(I)) * vely_(I); });
  });

  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcflux=[&](auto &u, auto f) {
      auto I_m = p.I - p.DI(2);
      auto I_p = p.I;
      auto lambda_m = 1.0;
      auto lambda_p = -1.0;
      auto var_m = u(I_m);
      auto var_p = u(I_p);
      auto flux_m = f(I_m);
      auto flux_p = f(I_p);
      return dAz * llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
    };

    fzrho_(p.I) = calcflux(rho_, [&](auto I) { return rho_(I) * velz_(I); });

    fzmomx_(p.I) = calcflux(momx_, [&](auto I) { return momx_(I) * velz_(I); });
    fzmomy_(p.I) = calcflux(momy_, [&](auto I) { return momy_(I) * velz_(I); });
    fzmomz_(p.I) = calcflux(
        momz_, [&](auto I) { return momz_(I) * velz_(I) + press_(I); });

    fzetot_(p.I) = calcflux(
        etot_, [&](auto I) { return (etot_(I) + press_(I)) * velz_(I); });
  });
}
#endif

#if 0
extern "C" void HydroToyCarpetX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Evolve;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dV1 = 1 / (dx * dy * dz);

  const Loop::GF3D1<const CCTK_REAL> rho_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                            rho_p);
  const Loop::GF3D1<const CCTK_REAL> momx_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momx_p);
  const Loop::GF3D1<const CCTK_REAL> momy_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momy_p);
  const Loop::GF3D1<const CCTK_REAL> momz_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momz_p);
  const Loop::GF3D1<const CCTK_REAL> etot_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             etot_p);

  const Loop::GF3D1<const CCTK_REAL> fxrho_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                            fxrho);
  const Loop::GF3D1<const CCTK_REAL> fxmomx_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                             fxmomx);
  const Loop::GF3D1<const CCTK_REAL> fxmomy_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                             fxmomy);
  const Loop::GF3D1<const CCTK_REAL> fxmomz_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                             fxmomz);
  const Loop::GF3D1<const CCTK_REAL> fxetot_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                             fxetot);

  const Loop::GF3D1<const CCTK_REAL> fyrho_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                            fyrho);
  const Loop::GF3D1<const CCTK_REAL> fymomx_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                             fymomx);
  const Loop::GF3D1<const CCTK_REAL> fymomy_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                             fymomy);
  const Loop::GF3D1<const CCTK_REAL> fymomz_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                             fymomz);
  const Loop::GF3D1<const CCTK_REAL> fyetot_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                             fyetot);

  const Loop::GF3D1<const CCTK_REAL> fzrho_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                            fzrho);
  const Loop::GF3D1<const CCTK_REAL> fzmomx_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                             fzmomx);
  const Loop::GF3D1<const CCTK_REAL> fzmomy_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                             fzmomy);
  const Loop::GF3D1<const CCTK_REAL> fzmomz_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                             fzmomz);
  const Loop::GF3D1<const CCTK_REAL> fzetot_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                             fzetot);

  const Loop::GF3D1<CCTK_REAL> rho_(cctkGH, {1, 1, 1}, {1, 1, 1}, rho);
  const Loop::GF3D1<CCTK_REAL> momx_(cctkGH, {1, 1, 1}, {1, 1, 1}, momx);
  const Loop::GF3D1<CCTK_REAL> momy_(cctkGH, {1, 1, 1}, {1, 1, 1}, momy);
  const Loop::GF3D1<CCTK_REAL> momz_(cctkGH, {1, 1, 1}, {1, 1, 1}, momz);
  const Loop::GF3D1<CCTK_REAL> etot_(cctkGH, {1, 1, 1}, {1, 1, 1}, etot);

  // Transport
  // dt rho + d_i (rho vel^i) = 0
  // dt mom_j + d_i (mom_j vel^i) = 0
  // dt etot + d_i (etot vel^i) = 0

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcupdate=[&](auto &fx, auto &fy, auto &fz) {
      return dV1 *
             ((fx(p.I + p.DI(0)) - fx(p.I)) + (fy(p.I + p.DI(1)) - fy(p.I)) +
              (fz(p.I + p.DI(2)) - fz(p.I)));
    };

    rho_(p.I) = rho_p_(p.I) - calcupdate(fxrho_, fyrho_, fzrho_);

    momx_(p.I) = momx_p_(p.I) - calcupdate(fxmomx_, fymomx_, fzmomx_);
    momy_(p.I) = momy_p_(p.I) - calcupdate(fxmomy_, fymomy_, fzmomy_);
    momz_(p.I) = momz_p_(p.I) - calcupdate(fxmomz_, fymomz_, fzmomz_);

    etot_(p.I) = etot_p_(p.I) - calcupdate(fxetot_, fyetot_, fzetot_);
  });
}
#endif

extern "C" void HydroToyCarpetX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D1<const CCTK_REAL> rho_(cctkGH, {1, 1, 1}, {1, 1, 1}, rho);
  const Loop::GF3D1<const CCTK_REAL> momx_(cctkGH, {1, 1, 1}, {1, 1, 1}, momx);
  const Loop::GF3D1<const CCTK_REAL> momy_(cctkGH, {1, 1, 1}, {1, 1, 1}, momy);
  const Loop::GF3D1<const CCTK_REAL> momz_(cctkGH, {1, 1, 1}, {1, 1, 1}, momz);
  const Loop::GF3D1<const CCTK_REAL> etot_(cctkGH, {1, 1, 1}, {1, 1, 1}, etot);

  const Loop::GF3D1<CCTK_REAL> regrid_error_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             regrid_error);

  if (false) {
    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      regrid_error_(p.I) = fabs(p.x) < 0.2;
    });
  }

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcerr = [&](auto &var_) {
      CCTK_REAL err{0};
      for (int d = 0; d < dim; ++d) {
        auto varm = var_(p.I - p.DI[d]);
        auto var0 = var_(p.I);
        auto varp = var_(p.I + p.DI[d]);
        err = fmax(err, fabs(varm - 2 * var0 + varp));
      }
      return err;
    };

    regrid_error_(p.I) = fmax5(calcerr(rho_), calcerr(momx_), calcerr(momy_),
                               calcerr(momz_), calcerr(etot_));
  });
}

// extern "C" void HydroToyCarpetX_Output(CCTK_ARGUMENTS) {
//  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Output;
//  DECLARE_CCTK_PARAMETERS;
//
//  const Loop::GF3D1<const CCTK_REAL> rho_(cctkGH, {1, 1, 1}, {1, 1, 1}, rho);
//  const Loop::GF3D1<const CCTK_REAL> momx_(cctkGH, {1, 1, 1}, {1, 1, 1},
//  momx); const Loop::GF3D1<const CCTK_REAL> momy_(cctkGH, {1, 1, 1}, {1, 1,
//  1}, momy); const Loop::GF3D1<const CCTK_REAL> momz_(cctkGH, {1, 1, 1}, {1,
//  1, 1}, momz); const Loop::GF3D1<const CCTK_REAL> etot_(cctkGH, {1, 1, 1},
//  {1, 1, 1}, etot);
//
//  const Loop::GF3D1<const CCTK_REAL> regrid_error_(cctkGH, {1, 1, 1}, {1, 1,
//  1},
//                                                   regrid_error);
//
//  const Loop::GF3D1<const CCTK_REAL> fxrho_(cctkGH, {0, 1, 1}, {0, 0, 0},
//                                            fxrho);
//  const Loop::GF3D1<const CCTK_REAL> fxmomx_(cctkGH, {0, 1, 1}, {0, 0, 0},
//                                             fxmomx);
//  const Loop::GF3D1<const CCTK_REAL> fxmomy_(cctkGH, {0, 1, 1}, {0, 0, 0},
//                                             fxmomy);
//  const Loop::GF3D1<const CCTK_REAL> fxmomz_(cctkGH, {0, 1, 1}, {0, 0, 0},
//                                             fxmomz);
//  const Loop::GF3D1<const CCTK_REAL> fxetot_(cctkGH, {0, 1, 1}, {0, 0, 0},
//                                             fxetot);
//
//  int levfac = cctk_levfac[0];
//  int lev = 0;
//  while (levfac > 1) {
//    levfac >>= 1;
//    lev += 1;
//  }
//
//#pragma omp critical(HydroToyCarpetX_Output)
//  {
//    const int oldprec = cout.precision(17);
//    cout << "iteration " << cctk_iteration << " level " << lev << ":\n";
//
//    Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
//      if (p.j == 0 && p.k == 0) {
//        cout << "[" << p.i << "," << p.j << "," << p.k << "] (" << p.x << ","
//             << p.y << "," << p.z << ") rho="
//             << rho_(p.I)
//             // << " momx=" << momx_(p.I) << " etot=" << etot_(p.I)
//             // << " err=" << regrid_error_(p.I)
//             << "\n";
//      }
//    });
//
//    Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
//      if (p.j == 0 && p.k == 0) {
//        cout << "[" << p.i << "," << p.j << "," << p.k << "] (" << p.x << ","
//             << p.y << "," << p.z << ") fxrho="
//             << fxrho_(p.I)
//             // << " fxmomx=" << fxmomx_(p.I) << " fxetot=" << fxetot_(p.I)
//             << "\n";
//      }
//    });
//
//    cout.precision(oldprec);
//  }
//}
//
} // namespace HydroToyCarpetX
