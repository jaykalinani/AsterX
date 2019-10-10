#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

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

  // const CCTK_REAL t = cctk_time;
  // const CCTK_REAL dt = CCTK_DELTA_TIME;

  if (CCTK_EQUALS(setup, "equilibrium")) {

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      rho[p.idx] = 1.0;
      momx[p.idx] = 0.0;
      momy[p.idx] = 0.0;
      momz[p.idx] = 0.0;
      etot[p.idx] = 1.0;
    });

  } else if (CCTK_EQUALS(setup, "sound wave")) {

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      rho[p.idx] = 1.0;
      momx[p.idx] = 0.0 + amplitude * sin(M_PI * p.x);
      momy[p.idx] = 0.0;
      momz[p.idx] = 0.0;
      etot[p.idx] = 1.0; // should add kinetic energy here
    });

  } else if (CCTK_EQUALS(setup, "shock tube")) {

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      if (p.x <= 0.0) {
        rho[p.idx] = 2.0;
        momx[p.idx] = 0.0;
        momy[p.idx] = 0.0;
        momz[p.idx] = 0.0;
        etot[p.idx] = 2.0;
      } else {
        rho[p.idx] = 1.0;
        momx[p.idx] = 0.0;
        momy[p.idx] = 0.0;
        momz[p.idx] = 0.0;
        etot[p.idx] = 1.0;
      }
    });

  } else if (CCTK_EQUALS(setup, "spherical shock")) {

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      CCTK_REAL r2 = pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2);
      if (r2 <= pow(shock_radius, 2)) {
        rho[p.idx] = 2.0;
        momx[p.idx] = 0.0;
        momy[p.idx] = 0.0;
        momz[p.idx] = 0.0;
        etot[p.idx] = 2.0;
      } else {
        rho[p.idx] = 1.0;
        momx[p.idx] = 0.0;
        momy[p.idx] = 0.0;
        momz[p.idx] = 0.0;
        etot[p.idx] = 1.0;
      }
    });

  } else {
    assert(0);
  }
}

extern "C" void HydroToyAMReX_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  // do nothing
}

extern "C" void HydroToyAMReX_Pressure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Equation of state: p = (gamma - 1) e

  // vel^j = delta^j_i mom_i / rho

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    CCTK_REAL ekin =
        sqrt(pow(momx[p.idx], 2) + pow(momy[p.idx], 2) + pow(momz[p.idx], 2)) /
        (2 * rho[p.idx]);
    CCTK_REAL eint = etot[p.idx] - ekin;
    press[p.idx] = (gamma - 1) * eint;

    velx[p.idx] = momx[p.idx] / rho[p.idx];
    vely[p.idx] = momy[p.idx] / rho[p.idx];
    velz[p.idx] = momz[p.idx] / rho[p.idx];
  });
}

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

extern "C" void HydroToyAMReX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> rho_(cctkGH, rho);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momx_(cctkGH, momx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momy_(cctkGH, momy);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momz_(cctkGH, momz);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> etot_(cctkGH, etot);

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> press_(cctkGH, press);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> velx_(cctkGH, velx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> vely_(cctkGH, vely);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> velz_(cctkGH, velz);

  const Loop::GF3D<CCTK_REAL, 0, 1, 1> fxrho_(cctkGH, fxrho);
  const Loop::GF3D<CCTK_REAL, 0, 1, 1> fxmomx_(cctkGH, fxmomx);
  const Loop::GF3D<CCTK_REAL, 0, 1, 1> fxmomy_(cctkGH, fxmomy);
  const Loop::GF3D<CCTK_REAL, 0, 1, 1> fxmomz_(cctkGH, fxmomz);
  const Loop::GF3D<CCTK_REAL, 0, 1, 1> fxetot_(cctkGH, fxetot);

  const Loop::GF3D<CCTK_REAL, 1, 0, 1> fyrho_(cctkGH, fyrho);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> fymomx_(cctkGH, fymomx);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> fymomy_(cctkGH, fymomy);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> fymomz_(cctkGH, fymomz);
  const Loop::GF3D<CCTK_REAL, 1, 0, 1> fyetot_(cctkGH, fyetot);

  const Loop::GF3D<CCTK_REAL, 1, 1, 0> fzrho_(cctkGH, fzrho);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> fzmomx_(cctkGH, fzmomx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> fzmomy_(cctkGH, fzmomy);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> fzmomz_(cctkGH, fzmomz);
  const Loop::GF3D<CCTK_REAL, 1, 1, 0> fzetot_(cctkGH, fzetot);

  // frho^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // fetot^i = (etot + press) vel^i

  Loop::loop_int<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcflux{[&](auto &u, auto f) {
      auto I_m = p.I - p.DI(0);
      auto I_p = p.I;
      auto lambda_m = 1.0;
      auto lambda_p = -1.0;
      auto var_m = u(I_m);
      auto var_p = u(I_p);
      auto flux_m = f(I_m);
      auto flux_p = f(I_p);
      return llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
    }};

    fxrho_(p.I) = calcflux(rho_, [&](auto I) { return rho_(I) * velx_(I); });

    fxmomx_(p.I) = calcflux(
        momx_, [&](auto I) { return momx_(I) * velx_(I) + press_(I); });
    fxmomy_(p.I) = calcflux(momy_, [&](auto I) { return momy_(I) * velx_(I); });
    fxmomz_(p.I) = calcflux(momz_, [&](auto I) { return momz_(I) * velx_(I); });

    fxetot_(p.I) = calcflux(
        etot_, [&](auto I) { return (etot_(I) + press_(I)) * velx_(I); });
  });

  Loop::loop_int<1, 0, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcflux{[&](auto &u, auto f) {
      auto I_m = p.I - p.DI(1);
      auto I_p = p.I;
      auto lambda_m = 1.0;
      auto lambda_p = -1.0;
      auto var_m = u(I_m);
      auto var_p = u(I_p);
      auto flux_m = f(I_m);
      auto flux_p = f(I_p);
      return llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
    }};

    fyrho_(p.I) = calcflux(rho_, [&](auto I) { return rho_(I) * vely_(I); });

    fymomx_(p.I) = calcflux(momx_, [&](auto I) { return momx_(I) * vely_(I); });
    fymomy_(p.I) = calcflux(
        momy_, [&](auto I) { return momy_(I) * vely_(I) + press_(I); });
    fymomz_(p.I) = calcflux(momz_, [&](auto I) { return momz_(I) * vely_(I); });

    fyetot_(p.I) = calcflux(
        etot_, [&](auto I) { return (etot_(I) + press_(I)) * vely_(I); });
  });

  Loop::loop_int<1, 1, 0>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcflux{[&](auto &u, auto f) {
      auto I_m = p.I - p.DI(2);
      auto I_p = p.I;
      auto lambda_m = 1.0;
      auto lambda_p = -1.0;
      auto var_m = u(I_m);
      auto var_p = u(I_p);
      auto flux_m = f(I_m);
      auto flux_p = f(I_p);
      return llf(lambda_m, lambda_p, var_m, var_p, flux_m, flux_p);
    }};

    fzrho_(p.I) = calcflux(rho_, [&](auto I) { return rho_(I) * velz_(I); });

    fzmomx_(p.I) = calcflux(momx_, [&](auto I) { return momx_(I) * velz_(I); });
    fzmomy_(p.I) = calcflux(momy_, [&](auto I) { return momy_(I) * velz_(I); });
    fzmomz_(p.I) = calcflux(
        momz_, [&](auto I) { return momz_(I) * velz_(I) + press_(I); });

    fzetot_(p.I) = calcflux(
        etot_, [&](auto I) { return (etot_(I) + press_(I)) * velz_(I); });
  });
}

extern "C" void HydroToyAMReX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> rho_p_(cctkGH, rho_p);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momx_p_(cctkGH, momx_p);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momy_p_(cctkGH, momy_p);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momz_p_(cctkGH, momz_p);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> etot_p_(cctkGH, etot_p);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxrho_(cctkGH, fxrho);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxmomx_(cctkGH, fxmomx);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxmomy_(cctkGH, fxmomy);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxmomz_(cctkGH, fxmomz);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxetot_(cctkGH, fxetot);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> fyrho_(cctkGH, fyrho);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> fymomx_(cctkGH, fymomx);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> fymomy_(cctkGH, fymomy);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> fymomz_(cctkGH, fymomz);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> fyetot_(cctkGH, fyetot);

  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> fzrho_(cctkGH, fzrho);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> fzmomx_(cctkGH, fzmomx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> fzmomy_(cctkGH, fzmomy);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> fzmomz_(cctkGH, fzmomz);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> fzetot_(cctkGH, fzetot);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> rho_(cctkGH, rho);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> momx_(cctkGH, momx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> momy_(cctkGH, momy);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> momz_(cctkGH, momz);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> etot_(cctkGH, etot);

  // Transport
  // dt rho + d_i (rho vel^i) = 0
  // dt mom_j + d_i (mom_j vel^i) = 0
  // dt etot + d_i (etot vel^i) = 0

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcupdate{[&](auto &fx, auto &fy, auto &fz) {
      return dt * ((fx(p.I + p.DI(0)) - fx(p.I)) / dx +
                   (fy(p.I + p.DI(1)) - fy(p.I)) / dy +
                   (fz(p.I + p.DI(2)) - fz(p.I)) / dz);
    }};

    rho_(p.I) = rho_p_(p.I) - calcupdate(fxrho_, fyrho_, fzrho_);

    momx_(p.I) = momx_p_(p.I) - calcupdate(fxmomx_, fymomx_, fzmomx_);
    momy_(p.I) = momy_p_(p.I) - calcupdate(fxmomy_, fymomy_, fzmomy_);
    momz_(p.I) = momz_p_(p.I) - calcupdate(fxmomz_, fymomz_, fzmomz_);

    etot_(p.I) = etot_p_(p.I) - calcupdate(fxetot_, fyetot_, fzetot_);
  });
}

extern "C" void HydroToyAMReX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> rho_(cctkGH, rho);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momx_(cctkGH, momx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momy_(cctkGH, momy);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momz_(cctkGH, momz);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> etot_(cctkGH, etot);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> regrid_error_(cctkGH, regrid_error);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    auto calcerr{[&](auto &var_) {
      CCTK_REAL err{0};
      for (int d = 0; d < dim; ++d) {
        auto varm = var_(p.I - p.DI(d));
        auto var0 = var_(p.I);
        auto varp = var_(p.I + p.DI(d));
        err = fmax(err, fabs(varm - 2 * var0 + varp));
      }
      return err;
    }};

    regrid_error_(p.I) =
        fmax(fmax(fmax(fmax(calcerr(rho_), calcerr(momx_)), calcerr(momy_)),
                  calcerr(momz_)),
             calcerr(etot_));
  });
}

extern "C" void HydroToyAMReX_Output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> rho_(cctkGH, rho);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momx_(cctkGH, momx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momy_(cctkGH, momy);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> momz_(cctkGH, momz);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> etot_(cctkGH, etot);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxrho_(cctkGH, fxrho);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxmomx_(cctkGH, fxmomx);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxmomy_(cctkGH, fxmomy);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxmomz_(cctkGH, fxmomz);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> fxetot_(cctkGH, fxetot);

  Loop::loop_all<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if (p.j == 0 && p.k == 0) {
#pragma omp critical(HydroToyAMReX_Output)
      {
        cout << "[" << p.i << "," << p.j << "," << p.k << "] (" << p.x << ","
             << p.y << "," << p.z << ") rho=" << rho_(p.I)
             << " momx=" << momx_(p.I) << " etot=" << etot_(p.I) << "\n";
      }
    }
  });

  Loop::loop_all<0, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if (p.j == 0 && p.k == 0) {
#pragma omp critical(HydroToyAMReX_Output)
      {
        cout << "[" << p.i << "," << p.j << "," << p.k << "] (" << p.x << ","
             << p.y << "," << p.z << ") fxrho=" << fxrho_(p.I)
             << " fxmomx=" << fxmomx_(p.I) << " fxetot=" << fxetot_(p.I)
             << "\n";
      }
    }
  });
}

} // namespace HydroToyAMReX
