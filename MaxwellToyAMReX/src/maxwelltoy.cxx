#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop.hxx>

#include <cassert>
#include <cmath>
#include <iostream>

namespace MaxwellToyAMReX {
using namespace std;

constexpr int dim = 3;

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct potential {
  T phi, ax, ay, az;
  potential operator-(const potential &y) const {
    const auto &x = *this;
    return {.phi = x.phi - y.phi,
            .ax = x.ax - y.ax,
            .ay = x.ay - y.ay,
            .az = x.az - y.az};
  }
  potential operator/(T a) const {
    const auto &x = *this;
    return {.phi = x.phi / a, .ax = x.ax / a, .ay = x.ay / a, .az = x.az / a};
  }
};

template <typename T> potential<T> plane_wave(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = 2 * M_PI * spatial_frequency_x;
  T ky = 2 * M_PI * spatial_frequency_y;
  T kz = 2 * M_PI * spatial_frequency_z;
  T omega = sqrt(pow(kx, 2) + pow(ky, 2) + pow(kz, 2));
  return {.phi = cos(omega * t - kx * x - ky * y - kz * z),
          .ax = omega / kx * cos(omega * t - kx * x - ky * y - kz * z),
          .ay = 0,
          .az = 0};
}

////////////////////////////////////////////////////////////////////////////////

// Derivatives
template <typename T, typename R> auto tderiv(R f(T t, T x, T y, T z), T dt) {
  return [=](T t, T x, T y, T z) {
    return (f(t + dt / 2, x, y, z) - f(t - dt / 2, x, y, z)) / dt;
  };
}
template <typename T, typename R> auto xderiv(R f(T t, T x, T y, T z), T dx) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x + dx / 2, y, z) - f(t, x - dx / 2, y, z)) / dx;
  };
}
template <typename T, typename R> auto yderiv(R f(T t, T x, T y, T z), T dy) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y + dy / 2, z) - f(t, x, y - dy / 2, z)) / dy;
  };
}
template <typename T, typename R> auto zderiv(R f(T t, T x, T y, T z), T dz) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z + dz / 2) - f(t, x, y, z - dz / 2)) / dz;
  };
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void MaxwellToyAMReX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL x0 = CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL y0 = CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL z0 = CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  if (CCTK_EQUALS(initial_condition, "plane wave")) {

    Loop::loop_all(cctkGH, [&](int i, int j, int k, int idx) {
      CCTK_REAL x = x0 + (cctk_lbnd[0] + i) * dx;
      CCTK_REAL y = y0 + (cctk_lbnd[1] + j) * dy;
      CCTK_REAL z = z0 + (cctk_lbnd[2] + k) * dz;

      phi[idx] = plane_wave(t + dt / 2, x, y, z).phi;
      ax[idx] = plane_wave(t, x + dx / 2, y, z).ax;
      ay[idx] = plane_wave(t, x, y + dy / 2, z).ay;
      az[idx] = plane_wave(t, x, y, z + dz / 2).az;

      ex[idx] =
          -xderiv(plane_wave<CCTK_REAL>, dx)(t + dt / 2, x + dx / 2, y, z).phi -
          tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, x + dx / 2, y, z).ax;
      ey[idx] =
          -yderiv(plane_wave<CCTK_REAL>, dy)(t + dt / 2, x, y + dy / 2, z).phi -
          tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, x, y + dy / 2, z).ay;
      ez[idx] =
          -zderiv(plane_wave<CCTK_REAL>, dz)(t + dt / 2, x, y, z + dz / 2).phi -
          tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, x, y, z + dz / 2).az;
      bx[idx] =
          yderiv(plane_wave<CCTK_REAL>, dy)(t, x, y + dy / 2, z + dz / 2).az -
          zderiv(plane_wave<CCTK_REAL>, dz)(t, x, y + dy / 2, z + dz / 2).ay;
      by[idx] =
          zderiv(plane_wave<CCTK_REAL>, dz)(t, x + dx / 2, y, z + dz / 2).ax -
          xderiv(plane_wave<CCTK_REAL>, dx)(t, x + dx / 2, y, z + dz / 2).az;
      bz[idx] =
          xderiv(plane_wave<CCTK_REAL>, dx)(t, x, y + dy / 2, z + dz / 2).ay -
          yderiv(plane_wave<CCTK_REAL>, dy)(t, x, y + dy / 2, z + dz / 2).ax;

      rho[idx] = 0.0;
      jx[idx] = 0.0;
      jy[idx] = 0.0;
      jz[idx] = 0.0;
    });

  } else {
    assert(0);
  }
}

extern "C" void MaxwellToyAMReX_Evolve1(CCTK_ARGUMENTS) {
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
    jx[idx] = jx_p[idx];
    jy[idx] = jy_p[idx];
    jz[idx] = jz_p[idx];

    bx[idx] = bx_p[idx] + dt * ((ey_p[idx + dk] - ey_p[idx]) / dz -
                                (ez_p[idx + dj] - ez_p[idx]) / dy);
    by[idx] = by_p[idx] + dt * ((ez_p[idx + di] - ez_p[idx]) / dx -
                                (ex_p[idx + dk] - ex_p[idx]) / dz);
    bz[idx] = bz_p[idx] + dt * ((ex_p[idx + dj] - ex_p[idx]) / dy -
                                (ey_p[idx + di] - ey_p[idx]) / dx);

    ax[idx] =
        ax_p[idx] - dt * ((phi_p[idx + di] - phi_p[idx]) / dx + ex_p[idx]);
    ay[idx] =
        ay_p[idx] - dt * ((phi_p[idx + dj] - phi_p[idx]) / dy + ey_p[idx]);
    az[idx] =
        az_p[idx] - dt * ((phi_p[idx + dk] - phi_p[idx]) / dz + ez_p[idx]);
  });
}

extern "C" void MaxwellToyAMReX_Evolve2(CCTK_ARGUMENTS) {
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
    rho[idx] = rho_p[idx] - dt * ((jx[idx] - jx[idx - di]) / dx +
                                  (jy[idx] - jy[idx - dj]) / dy +
                                  (jz[idx] - jz[idx - dk]) / dz);

    ex[idx] = ex_p[idx] - dt * ((by[idx] - by[idx - dk]) / dz -
                                (bz[idx] - bz[idx - dj]) / dy + jx[idx]);
    ey[idx] = ey_p[idx] - dt * ((bz[idx] - bz[idx - di]) / dx -
                                (bx[idx] - bx[idx - dk]) / dz + jy[idx]);
    ez[idx] = ez_p[idx] - dt * ((bx[idx] - bx[idx - dj]) / dy -
                                (by[idx] - by[idx - di]) / dx + jz[idx]);

    phi[idx] = phi_p[idx] - dt * ((ax[idx] - ax[idx - di]) / dx +
                                  (ay[idx] - ay[idx - dj]) / dy +
                                  (ay[idx] - az[idx - dk]) / dz);
  });
}

extern "C" void MaxwellToyAMReX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  Loop::loop_int(cctkGH,
                 [&](int i, int j, int k, int idx) { regrid_error[idx] = 0; });
}

extern "C" void MaxwellToyAMReX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL x0 = CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL y0 = CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL z0 = CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  if (CCTK_EQUALS(initial_condition, "plane wave")) {

    Loop::loop_all(cctkGH, [&](int i, int j, int k, int idx) {
      CCTK_REAL x = x0 + (cctk_lbnd[0] + i) * dx;
      CCTK_REAL y = y0 + (cctk_lbnd[1] + j) * dy;
      CCTK_REAL z = z0 + (cctk_lbnd[2] + k) * dz;

      phierr[idx] = phi[idx] - plane_wave(t + dt / 2, x, y, z).phi;
      axerr[idx] = ax[idx] - plane_wave(t, x + dx / 2, y, z).ax;
      ayerr[idx] = ay[idx] - plane_wave(t, x, y + dy / 2, z).ay;
      azerr[idx] = az[idx] - plane_wave(t, x, y, z + dz / 2).az;

      exerr[idx] =
          ex[idx] -
          (-xderiv(plane_wave<CCTK_REAL>, dx)(t + dt / 2, x + dx / 2, y, z)
                .phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, x + dx / 2, y, z).ax);
      eyerr[idx] =
          ey[idx] -
          (-yderiv(plane_wave<CCTK_REAL>, dy)(t + dt / 2, x, y + dy / 2, z)
                .phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, x, y + dy / 2, z).ay);
      ezerr[idx] =
          ez[idx] -
          (-zderiv(plane_wave<CCTK_REAL>, dz)(t + dt / 2, x, y, z + dz / 2)
                .phi -
           tderiv(plane_wave<CCTK_REAL>, dt)(t + dt / 2, x, y, z + dz / 2).az);
      bxerr[idx] =
          bx[idx] -
          (yderiv(plane_wave<CCTK_REAL>, dy)(t, x, y + dy / 2, z + dz / 2).az -
           zderiv(plane_wave<CCTK_REAL>, dz)(t, x, y + dy / 2, z + dz / 2).ay);
      byerr[idx] =
          by[idx] -
          (zderiv(plane_wave<CCTK_REAL>, dz)(t, x + dx / 2, y, z + dz / 2).ax -
           xderiv(plane_wave<CCTK_REAL>, dx)(t, x + dx / 2, y, z + dz / 2).az);
      bzerr[idx] =
          bz[idx] -
          (xderiv(plane_wave<CCTK_REAL>, dx)(t, x, y + dy / 2, z + dz / 2).ay -
           yderiv(plane_wave<CCTK_REAL>, dy)(t, x, y + dy / 2, z + dz / 2).ax);

      rhoerr[idx] = 0.0;
      jxerr[idx] = 0.0;
      jyerr[idx] = 0.0;
      jzerr[idx] = 0.0;
    });

  } else {
    assert(0);
  }
}

} // namespace MaxwellToyAMReX
