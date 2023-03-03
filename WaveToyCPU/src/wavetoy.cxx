#include <fixmath.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop.hxx>

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace WaveToyCPU {
using namespace std;
using namespace Loop;

// Linear interpolation between (i0, x0) and (i1, x1)
template <typename Y, typename X> Y linterp(Y y0, Y y1, X x0, X x1, X x) {
  return Y(x - x0) / Y(x1 - x0) * y0 + Y(x - x1) / Y(x0 - x1) * y1;
}

// Spline with compact support of radius 1 and volume 1
template <typename T> T spline(T r) {
  if (r >= 1.0)
    return 0.0;
  constexpr CCTK_REAL f = dim == 1   ? 1.0
                          : dim == 2 ? 24.0 / 7.0 / M_PI
                          : dim == 3 ? 4.0 / M_PI
                                     : -1;
  const T r2 = pow(r, 2);
  return f * (r <= 0.5 ? 1 - 2 * r2 : 2 + r * (-4 + 2 * r));
}

// The potential for the spline
template <typename T> T spline_potential(T r) {
  // \Laplace u = 4 \pi \rho
  if (r >= 1.0)
    return -1 / r;
  static_assert(dim == 3, "");
  const T r2 = pow(r, 2);
  return r <= 0.5
             ? -7 / T(3) + r2 * (8 / T(3) - 8 / T(5) * r2)
             : (1 / T(15) +
                r * (-40 / T(15) +
                     r2 * (80 / T(15) + r * (-80 / T(15) + r * 24 / T(15))))) /
                   r;
}

// Time derivative
// TODO: Use dual numbers for derivative
template <typename F, typename T> auto timederiv(const F &f, T dt) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z) - f(t - dt, x, y, z)) / dt;
  };
}

// Gradient
template <typename T> auto xderiv(T f(T t, T x, T y, T z), T dx) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x + dx, y, z) - f(t, x - dx, y, z)) / (2 * dx);
  };
}
template <typename T> auto yderiv(T f(T t, T x, T y, T z), T dy) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y + dy, z) - f(t, x, y - dy, z)) / (2 * dy);
  };
}
template <typename T> auto zderiv(T f(T t, T x, T y, T z), T dz) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z + dz) - f(t, x, y, z - dz)) / (2 * dz);
  };
}

////////////////////////////////////////////////////////////////////////////////

// Standing wave
CCTK_REAL standing(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto kx = 2 * M_PI * spatial_frequency_x;
  auto ky = 2 * M_PI * spatial_frequency_y;
  auto kz = 2 * M_PI * spatial_frequency_z;
  auto omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
}

// Periodic Gaussian
CCTK_REAL periodic_gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                            CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto kx = M_PI * spatial_frequency_x;
  auto ky = M_PI * spatial_frequency_y;
  auto kz = M_PI * spatial_frequency_z;
  auto omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return exp(-0.5 * pow(sin(kx * x + ky * y + kz * z - omega * t) / width, 2));
}

// Gaussian
CCTK_REAL gaussian(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  // u(t,r) = (f(r-t) - f(r+t)) / r
  auto r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  auto f = [&](auto x) { return exp(-0.5 * pow(x / width, 2)); };
  auto fx = [&](auto x) { return -x / pow(width, 2) * f(x); };
  if (r < 1.0e-8)
    // Use L'Hôpital's rule for small r
    return fx(r - t) - fx(r + t);
  else
    return (f(r - t) - f(r + t)) / r;
}

// Central potential
CCTK_REAL central_potential(CCTK_REAL t, CCTK_REAL x, CCTK_REAL y,
                            CCTK_REAL z) {
  DECLARE_CCTK_PARAMETERS;
  auto r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
  return -central_point_charge / pow(central_point_radius, dim) *
         spline_potential(r / central_point_radius);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyCPU_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<CCTK_REAL> gf_psi(layout, psi);

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = standing(t, p.x, p.y, p.z);
      gf_psi(p.I) = timederiv(standing, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = periodic_gaussian(t, p.x, p.y, p.z);
      gf_psi(p.I) = timederiv(periodic_gaussian, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = gaussian(t, p.x, p.y, p.z);
      gf_psi(p.I) = timederiv(gaussian, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "central potential")) {

    loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phi(p.I) = central_potential(t, p.x, p.y, p.z);
      gf_psi(p.I) = timederiv(central_potential, dt)(t, p.x, p.y, p.z);
    });

  } else {
    assert(0);
  }
}

extern "C" void WaveToyCPU_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Sync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyCPU_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<const CCTK_REAL> gf_psi(layout, psi);
  const GF3D2<CCTK_REAL> gf_regrid_error(layout, regrid_error);

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    CCTK_REAL base_phi = fabs(gf_phi(p.I)) + fabs(phi_abs);
    CCTK_REAL errx_phi =
        fabs(gf_phi(p.I - p.DI[0]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[0])) /
        base_phi;
    CCTK_REAL erry_phi =
        fabs(gf_phi(p.I - p.DI[1]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[1])) /
        base_phi;
    CCTK_REAL errz_phi =
        fabs(gf_phi(p.I - p.DI[2]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[2])) /
        base_phi;
    CCTK_REAL base_psi = fabs(gf_psi(p.I)) + fabs(psi_abs);
    CCTK_REAL errx_psi =
        fabs(gf_psi(p.I - p.DI[0]) - 2 * gf_psi(p.I) + gf_psi(p.I + p.DI[0])) /
        base_psi;
    CCTK_REAL erry_psi =
        fabs(gf_psi(p.I - p.DI[1]) - 2 * gf_psi(p.I) + gf_psi(p.I + p.DI[1])) /
        base_psi;
    CCTK_REAL errz_psi =
        fabs(gf_psi(p.I - p.DI[2]) - 2 * gf_psi(p.I) + gf_psi(p.I + p.DI[2])) /
        base_psi;
    gf_regrid_error(p.I) =
        errx_phi + erry_phi + errz_phi + errx_psi + erry_psi + errz_psi;
  });
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyCPU_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_RHS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<const CCTK_REAL> gf_psi(layout, psi);
  const GF3D2<CCTK_REAL> gf_phirhs(layout, phirhs);
  const GF3D2<CCTK_REAL> gf_psirhs(layout, psirhs);

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    CCTK_REAL ddx_phi =
        (gf_phi(p.I - p.DI[0]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[0])) /
        pow(p.dx, 2);
    CCTK_REAL ddy_phi =
        (gf_phi(p.I - p.DI[1]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[1])) /
        pow(p.dy, 2);
    CCTK_REAL ddz_phi =
        (gf_phi(p.I - p.DI[2]) - 2 * gf_phi(p.I) + gf_phi(p.I + p.DI[2])) /
        pow(p.dz, 2);
    gf_psirhs(p.I) = ddx_phi + ddy_phi + ddz_phi - pow(mass, 2) * gf_phi(p.I) +
                     4 * M_PI * central_potential(t, p.x, p.y, p.z);
    gf_phirhs(p.I) = gf_psi(p.I);
  });
}

extern "C" void WaveToyCPU_RHSSync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_RHSSync;
  DECLARE_CCTK_PARAMETERS;

  // Do nothing
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyCPU_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Energy;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<const CCTK_REAL> gf_psi(layout, psi);
  const GF3D2<CCTK_REAL> gf_eps(layout, eps);

  loop_int<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
    CCTK_REAL dt_phi = gf_psi(p.I);
    CCTK_REAL dx_phi =
        (gf_phi(p.I + p.DI[0]) - gf_phi(p.I - p.DI[0])) / (2 * p.dx);
    CCTK_REAL dy_phi =
        (gf_phi(p.I + p.DI[1]) - gf_phi(p.I - p.DI[1])) / (2 * p.dy);
    CCTK_REAL dz_phi =
        (gf_phi(p.I + p.DI[2]) - gf_phi(p.I - p.DI[2])) / (2 * p.dz);
    gf_eps(p.I) =
        (pow(dt_phi, 2) + pow(dx_phi, 2) + pow(dy_phi, 2) + pow(dz_phi, 2)) / 2;
  });
}

extern "C" void WaveToyCPU_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_WaveToyCPU_Error;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;

  const array<int, dim> indextype = {1, 1, 1};
  const GF3D2layout layout(cctkGH, indextype);
  const GF3D2<const CCTK_REAL> gf_phi(layout, phi);
  const GF3D2<const CCTK_REAL> gf_psi(layout, psi);
  const GF3D2<CCTK_REAL> gf_phierr(layout, phierr);
  const GF3D2<CCTK_REAL> gf_psierr(layout, psierr);

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phierr(p.I) = gf_phi(p.I) - standing(t, p.x, p.y, p.z);
      gf_psierr(p.I) = gf_psi(p.I) - timederiv(standing, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phierr(p.I) = gf_phi(p.I) - periodic_gaussian(t, p.x, p.y, p.z);
      gf_psierr(p.I) =
          gf_psi(p.I) - timederiv(periodic_gaussian, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "Gaussian")) {

    loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phierr(p.I) = gf_phi(p.I) - gaussian(t, p.x, p.y, p.z);
      gf_psierr(p.I) = gf_psi(p.I) - timederiv(gaussian, dt)(t, p.x, p.y, p.z);
    });

  } else if (CCTK_EQUALS(initial_condition, "central potential")) {

    loop_all<1, 1, 1>(cctkGH, [&](const PointDesc &p) {
      gf_phierr(p.I) = gf_phi(p.I) - central_potential(t, p.x, p.y, p.z);
      gf_psierr(p.I) =
          gf_psi(p.I) - timederiv(central_potential, dt)(t, p.x, p.y, p.z);
    });

  } else {
    assert(0);
  }
}

} // namespace WaveToyCPU
