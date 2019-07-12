#include <AMReX.hxx>
using namespace AMReX;

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX_PlotFileUtil.H>
using namespace amrex;

#include <omp.h>

#include <cmath>
using namespace std;

namespace WaveToyAMReX {

constexpr int dim = 3;

// Linear interpolation between (i0, x0) and (i1, x1)
template <typename Y, typename X> Y linterp(Y y0, Y y1, X x0, X x1, X x) {
  return Y(x - x0) / Y(x1 - x0) * y0 + Y(x - x1) / Y(x0 - x1) * y1;
}

// Square
template <typename T> T sqr(T x) { return x * x; }

// Standing wave
template <typename T> T standing(T t, T x, T y, T z) {
  T Lx = 2, Ly = 2, Lz = 2;
  T kx = 2 * M_PI / Lx, ky = 2 * M_PI / Ly, kz = 2 * M_PI / Lz;
  T omega = sqrt(sqr(kx) + sqr(ky) + sqr(kz));
  return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void WaveToyAMReX_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL *restrict const x0 = ghext->geom.ProbLo();
  const CCTK_REAL *restrict const x1 = ghext->geom.ProbHi();
  for (int d = 0; d < dim; ++d) {
    CCTK_REAL dx = (x1[d] - x0[d]) / ghext->ncells;
    cctkGH->cctk_origin_space[d] = x0[d] + 0.5 * dx;
    cctkGH->cctk_delta_space[d] = dx;
  }
  CCTK_VINFO("x0=[%g,%g,%g]", cctkGH->cctk_origin_space[0],
             cctkGH->cctk_origin_space[1], cctkGH->cctk_origin_space[2]);
  CCTK_VINFO("dx=[%g,%g,%g]", cctkGH->cctk_delta_space[0],
             cctkGH->cctk_delta_space[1], cctkGH->cctk_delta_space[2]);

  CCTK_REAL mindx = 1.0 / 0.0;
  for (int d = 0; d < dim; ++d)
    mindx = fmin(mindx, cctkGH->cctk_delta_space[d]);

  cctkGH->cctk_time = 0.0;
  cctkGH->cctk_delta_time = 0.5 * mindx;
  CCTK_VINFO("t=%g", cctkGH->cctk_time);
  CCTK_VINFO("dt=%g", cctkGH->cctk_delta_time);
}

extern "C" void WaveToyAMReX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL dt = cctk_delta_time;
  const CCTK_REAL *restrict const x0 = cctk_origin_space;
  const CCTK_REAL *restrict const dx = cctk_delta_space;

  for (int k = 0; k < cctk_lsh[2]; ++k)
    for (int j = 0; j < cctk_lsh[1]; ++j)
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
        CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
        CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
        phi[idx] = standing(t, x, y, z);
        phi_p[idx] = standing(t - dt, x, y, z);
      }
}

extern "C" void WaveToyAMReX_Cycle(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Cycle time levels
  MultiFab::Copy(ghext->mfab, ghext->mfab, 1, 2, 1, ghext->nghostzones);
  MultiFab::Copy(ghext->mfab, ghext->mfab, 0, 1, 1, ghext->nghostzones);

  // Step time
  // cctkGH->cctk_time += cctkGH->cctk_delta_time;
  CCTK_VINFO("t=%g", cctkGH->cctk_time);
}

extern "C" void WaveToyAMReX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL *restrict const dx = cctk_delta_space;
  const CCTK_REAL dt = cctk_delta_time;

  constexpr int di = 1;
  const int dj = di * cctk_ash[0];
  const int dk = dj * cctk_ash[1];

  for (int k = 0; k < cctk_lsh[2]; ++k)
    for (int j = 0; j < cctk_lsh[1]; ++j)
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL ddx_phi =
            (phi_p[idx - di] - 2 * phi_p[idx] + phi_p[idx + di]) / sqr(dx[0]);
        CCTK_REAL ddy_phi =
            (phi_p[idx - dj] - 2 * phi_p[idx] + phi_p[idx + dj]) / sqr(dx[1]);
        CCTK_REAL ddz_phi =
            (phi_p[idx - dk] - 2 * phi_p[idx] + phi_p[idx + dk]) / sqr(dx[2]);
        phi[idx] = -phi_p_p[idx] + 2 * phi_p[idx] +
                   sqr(dt) * (ddx_phi + ddy_phi + ddz_phi);
      }
}

extern "C" void WaveToyAMReX_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Synchronize
  ghext->mfab.FillBoundary(ghext->geom.periodicity());
}

extern "C" void WaveToyAMReX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL t = cctk_time;
  const CCTK_REAL *restrict const x0 = cctk_origin_space;
  const CCTK_REAL *restrict const dx = cctk_delta_space;

  for (int k = 0; k < cctk_lsh[2]; ++k)
    for (int j = 0; j < cctk_lsh[1]; ++j)
#pragma omp simd
      for (int i = 0; i < cctk_lsh[0]; ++i) {
        const int idx = CCTK_GFINDEX3D(cctkGH, i, j, k);
        CCTK_REAL x = x0[0] + (cctk_lbnd[0] + i) * dx[0];
        CCTK_REAL y = x0[1] + (cctk_lbnd[1] + j) * dx[1];
        CCTK_REAL z = x0[2] + (cctk_lbnd[2] + k) * dx[2];
        err[idx] = phi[idx] - standing(t, x, y, z);
      }
}

extern "C" void WaveToyAMReX_Output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Output phi
  string filename = amrex::Concatenate("wavetoy/phi", cctk_iteration, 6);
  WriteSingleLevelPlotfile(filename, ghext->mfab,
                           {"phi", "phi_p", "phi_p_p", "err"}, ghext->geom,
                           cctk_time, cctk_iteration);
}

} // namespace WaveToyAMReX
