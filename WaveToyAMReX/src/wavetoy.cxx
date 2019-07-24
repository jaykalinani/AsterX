#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

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
    // return (f(t + dt, x, y, z) - f(t, x, y, z)) / dt;
    return (f(t, x, y, z) - f(t - dt, x, y, z)) / dt;
  };
}

// Standing wave
template <typename T> T standing(T t, T x, T y, T z) {
  T kx = 2 * M_PI / 2, ky = 2 * M_PI / 2, kz = 2 * M_PI / 2;
  // T kx = 2 * M_PI / 2, ky = 0, kz = 0;
  T omega = sqrt(sqr(kx) + sqr(ky) + sqr(kz));
  return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
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

  for (int k = imin[2]; k < imax[2]; ++k) {
    for (int j = imin[1]; j < imax[1]; ++j) {
#pragma omp simd
      for (int i = imin[0]; i < imax[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
        CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
        CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
        phi[idx] = standing(t, x, y, z);
        phi_p[idx] = standing(t - dt, x, y, z);
        // psi[idx] = timederiv(standing, dt)(t, x, y, z);
      }
    }
  }
}

extern "C" void WaveToyAMReX_Tag(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL x0[dim], dx[dim];
  for (int d = 0; d < dim; ++d) {
    x0[d] = CCTK_ORIGIN_SPACE(d);
    dx[d] = CCTK_DELTA_SPACE(d);
  }

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
        CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
        CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
        CCTK_REAL r = sqrt(sqr(x) + sqr(y) + sqr(z));
        // CCTK_REAL r = fabs(x);
        tag[idx] = r <= 0.5 / cctk_levfac[0] - 0.5 * dx[0];
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
        // CCTK_REAL ddx_psi =
        //     (psi_p[idx - di] - 2 * psi_p[idx] + psi_p[idx + di]) /
        //     sqr(dx[0]);
        // CCTK_REAL ddy_psi =
        //     (psi_p[idx - dj] - 2 * psi_p[idx] + psi_p[idx + dj]) /
        //     sqr(dx[1]);
        // CCTK_REAL ddz_psi =
        //     (psi_p[idx - dk] - 2 * psi_p[idx] + psi_p[idx + dk]) /
        //     sqr(dx[2]);
        // phi[idx] = phi_p[idx] + dt * psi_p[idx];
        // psi[idx] = psi_p[idx] + dt * (ddx_phi + ddy_phi + ddz_phi) +
        //            sqr(dt) * (ddx_psi + ddy_psi + ddz_psi);
        // phi[idx] = phi_p[idx] + dt * psi_p[idx];
        // psi[idx] = psi_p[idx] + dt * (ddx_phi + ddy_phi + ddz_phi);
        phi[idx] = 2 * phi_p[idx] - phi_p_p[idx] +
                   sqr(dt) * (ddx_phi + ddy_phi + ddz_phi);
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
        CCTK_REAL dt_phi = (phi[idx] - phi_p[idx]) / dt;
        CCTK_REAL dx_phi = (phi[idx + di] - phi[idx - di]) / (2 * dx[0]);
        CCTK_REAL dy_phi = (phi[idx + dj] - phi[idx - dj]) / (2 * dx[1]);
        CCTK_REAL dz_phi = (phi[idx + dk] - phi[idx - dk]) / (2 * dx[2]);
        eps[idx] = (sqr(dt_phi) + sqr(dx_phi) + sqr(dy_phi) + sqr(dz_phi)) / 2;
      }
    }
  }
}

extern "C" void WaveToyAMReX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  // const CCTK_REAL dt = CCTK_DELTA_TIME;
  CCTK_REAL x0[dim], dx[dim];
  for (int d = 0; d < dim; ++d) {
    x0[d] = CCTK_ORIGIN_SPACE(d);
    dx[d] = CCTK_DELTA_SPACE(d);
  }

  for (int k = 0; k < cctk_lsh[2]; ++k) {
    for (int j = 0; j < cctk_lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
        CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
        CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
        phierr[idx] = phi[idx] - standing(t, x, y, z);
        // psierr[idx] = psi[idx] - timederiv(standing, dt)(t, x, y, z);
      }
    }
  }
}

} // namespace WaveToyAMReX
