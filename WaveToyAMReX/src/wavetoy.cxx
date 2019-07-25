#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <iostream>

namespace WaveToyAMReX {
using namespace std;

constexpr int dim = 3;

// Linear interpolation between (i0, x0) and (i1, x1)
template <typename Y, typename X> Y linterp(Y y0, Y y1, X x0, X x1, X x) {
  return Y(x - x0) / Y(x1 - x0) * y0 + Y(x - x1) / Y(x0 - x1) * y1;
}

// Square
template <typename T> T sqr(T x) { return x * x; }

// Time derivative
template <typename T> auto timederiv(T f(T t, T x, T y, T z), T dt) {
  return [=](T t, T x, T y, T z) {
    return (f(t, x, y, z) - f(t - dt, x, y, z)) / dt;
  };
}

// Standing wave
template <typename T> T standing(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = 2 * M_PI * spatial_frequency_x;
  T ky = 2 * M_PI * spatial_frequency_y;
  T kz = 2 * M_PI * spatial_frequency_z;
  T omega = sqrt(sqr(kx) + sqr(ky) + sqr(kz));
  return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
}

// Standing wave
template <typename T> T gaussian(T t, T x, T y, T z) {
  DECLARE_CCTK_PARAMETERS;
  T kx = M_PI * spatial_frequency_x;
  T ky = M_PI * spatial_frequency_y;
  T kz = M_PI * spatial_frequency_z;
  T omega = sqrt(sqr(kx) + sqr(ky) + sqr(kz));
  return exp(-0.5 * sqr(sin(kx * x + ky * y + kz * z - omega * t) / width));
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyAMReX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  CCTK_REAL x0[dim], dx[dim];
  for (int d = 0; d < dim; ++d) {
    x0[d] = CCTK_ORIGIN_SPACE(d);
    dx[d] = CCTK_DELTA_SPACE(d);
  }

  int imin[dim], imax[dim];
  for (int d = 0; d < dim; ++d) {
    imin[d] = cctk_bbox[2 * d] ? 0 : cctk_nghostzones[d];
    imax[d] = cctk_lsh[d] - (cctk_bbox[2 * d + 1] ? 0 : cctk_nghostzones[d]);
  }

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
        for (int i = imin[0]; i < imax[0]; ++i) {
          const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
          CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
          CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
          CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
          phi[idx] = standing(t, x, y, z);
          // phi_p[idx] = standing(t - dt, x, y, z);
          psi[idx] = timederiv(standing, dt)(t, x, y, z);
        }
      }
    }

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
        for (int i = imin[0]; i < imax[0]; ++i) {
          const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
          CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
          CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
          CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
          phi[idx] = gaussian(t, x, y, z);
          // phi_p[idx] = gaussian(t - dt, x, y, z);
          psi[idx] = timederiv(gaussian, dt)(t, x, y, z);
        }
      }
    }

  } else {
    assert(0);
  }
}

extern "C" void WaveToyAMReX_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  constexpr int di = 1;
  const int dj = di * cctk_ash[0];
  const int dk = dj * cctk_ash[1];

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL base_phi = fabs(phi[idx]) + fabs(phi_abs);
        CCTK_REAL errx_phi =
            fabs(phi[idx - di] - 2 * phi[idx] + phi[idx + di]) / base_phi;
        CCTK_REAL erry_phi =
            fabs(phi[idx - dj] - 2 * phi[idx] + phi[idx + dj]) / base_phi;
        CCTK_REAL errz_phi =
            fabs(phi[idx - dk] - 2 * phi[idx] + phi[idx + dk]) / base_phi;
        regrid_error[idx] = errx_phi + erry_phi + errz_phi;
      }
    }
  }
}

extern "C" void WaveToyAMReX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  CCTK_REAL dx[dim];
  for (int d = 0; d < dim; ++d) {
    dx[d] = CCTK_DELTA_SPACE(d);
  }

  int imin[dim], imax[dim];
  for (int d = 0; d < dim; ++d) {
    imin[d] = cctk_bbox[2 * d] ? 0 : cctk_nghostzones[d];
    imax[d] = cctk_lsh[d] - (cctk_bbox[2 * d + 1] ? 0 : cctk_nghostzones[d]);
  }

  constexpr int di = 1;
  const int dj = di * cctk_ash[0];
  const int dk = dj * cctk_ash[1];

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
      for (int i = imin[0]; i < imax[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL ddx_phi =
            (phi_p[idx - di] - 2 * phi_p[idx] + phi_p[idx + di]) / sqr(dx[0]);
        CCTK_REAL ddy_phi =
            (phi_p[idx - dj] - 2 * phi_p[idx] + phi_p[idx + dj]) / sqr(dx[1]);
        CCTK_REAL ddz_phi =
            (phi_p[idx - dk] - 2 * phi_p[idx] + phi_p[idx + dk]) / sqr(dx[2]);
        // phi[idx] = 2 * phi_p[idx] - phi_p_p[idx] +
        //            sqr(dt) * (ddx_phi + ddy_phi + ddz_phi);
        psi[idx] = psi_p[idx] + dt * (ddx_phi + ddy_phi + ddz_phi);
        phi[idx] = phi_p[idx] + dt * psi[idx];
      }
    }
  }
}

extern "C" void WaveToyAMReX_Energy(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  CCTK_REAL dx[dim];
  for (int d = 0; d < dim; ++d)
    dx[d] = CCTK_DELTA_SPACE(d);

  int imin[dim], imax[dim];
  for (int d = 0; d < dim; ++d) {
    imin[d] = cctk_bbox[2 * d] ? 0 : cctk_nghostzones[d];
    imax[d] = cctk_lsh[d] - (cctk_bbox[2 * d + 1] ? 0 : cctk_nghostzones[d]);
  }

  constexpr int di = 1;
  const int dj = di * cctk_ash[0];
  const int dk = dj * cctk_ash[1];

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
      for (int i = imin[0]; i < imax[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL ddx_phi =
            (phi[idx - di] - 2 * phi[idx] + phi[idx + di]) / sqr(dx[0]);
        CCTK_REAL ddy_phi =
            (phi[idx - dj] - 2 * phi[idx] + phi[idx + dj]) / sqr(dx[1]);
        CCTK_REAL ddz_phi =
            (phi[idx - dk] - 2 * phi[idx] + phi[idx + dk]) / sqr(dx[2]);
        CCTK_REAL psi_n = psi[idx] + dt * (ddx_phi + ddy_phi + ddz_phi);
        CCTK_REAL dt_phi_p = psi_n;
        CCTK_REAL dt_phi_m = psi_p[idx];
        CCTK_REAL dx_phi_p = (phi[idx + di] - phi[idx]) / dx[0];
        CCTK_REAL dx_phi_m = (phi[idx] - phi[idx - di]) / dx[0];
        CCTK_REAL dy_phi_p = (phi[idx + dj] - phi[idx]) / dx[1];
        CCTK_REAL dy_phi_m = (phi[idx] - phi[idx - dj]) / dx[1];
        CCTK_REAL dz_phi_p = (phi[idx + dk] - phi[idx]) / dx[2];
        CCTK_REAL dz_phi_m = (phi[idx] - phi[idx - dk]) / dx[2];
        eps[idx] =
            (sqr(dt_phi_p) + sqr(dt_phi_m) + sqr(dx_phi_p) + sqr(dx_phi_m) +
             sqr(dy_phi_p) + sqr(dy_phi_m) + sqr(dz_phi_p) + sqr(dz_phi_m)) /
            4;
      }
    }
  }
}

extern "C" void WaveToyAMReX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = CCTK_DELTA_TIME;
  CCTK_REAL x0[dim], dx[dim];
  for (int d = 0; d < dim; ++d) {
    x0[d] = CCTK_ORIGIN_SPACE(d);
    dx[d] = CCTK_DELTA_SPACE(d);
  }

  if (CCTK_EQUALS(initial_condition, "standing wave")) {

    for (int k = 0; k < cctk_lsh[2]; ++k) {
      for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
        for (int i = 0; i < cctk_lsh[0]; ++i) {
          const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
          CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
          CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
          CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
          phierr[idx] = phi[idx] - standing(t, x, y, z);
          psierr[idx] = psi[idx] - timederiv(standing, dt)(t, x, y, z);
        }
      }
    }

  } else if (CCTK_EQUALS(initial_condition, "periodic Gaussian")) {

    for (int k = 0; k < cctk_lsh[2]; ++k) {
      for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
        for (int i = 0; i < cctk_lsh[0]; ++i) {
          const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
          CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
          CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
          CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
          phierr[idx] = phi[idx] - gaussian(t, x, y, z);
          psierr[idx] = psi[idx] - timederiv(gaussian, dt)(t, x, y, z);
        }
      }
    }
  }
}

} // namespace WaveToyAMReX
