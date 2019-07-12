#include <AMReX.hxx>
using namespace AMReX;

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX_PlotFileUtil.H>
using namespace amrex;

#include <atomic>
#include <cmath>
using namespace std;

namespace WaveToyAMReX {

// Linear interpolation between (i0, x0) and (i1, x1)
template <typename Y, typename X> Y linterp(Y y0, Y y1, X x0, X x1, X x) {
  return Y(x - x0) / Y(x1 - x0) * y0 + Y(x - x1) / Y(x0 - x1) * y1;
}

// Square
template <typename T> T sqr(T x) { return x * x; }

// Standing wave
template <typename T> T standing(T t, T x, T y, T z) {
  T Lx = 4, Ly = 4, Lz = 4;
  T kx = 2 * M_PI / Lx, ky = 2 * M_PI / Ly, kz = 2 * M_PI / Lz;
  T omega = sqrt(sqr(kx) + sqr(ky) + sqr(kz));
  return cos(omega * t) * cos(kx * x) * cos(ky * y) * cos(kz * z);
}

extern "C" void WaveToyAMReX_Check(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_VINFO("Check iterators");

  {
    int count = 0;
    for (MFIter mfi(ghext->mfab); mfi.isValid(); ++mfi) {
      const Box &fbx = mfi.fabbox();
      const Box &bx = mfi.validbox();

      const Dim3 amin = lbound(fbx);
      const Dim3 amax = ubound(fbx);
      const Dim3 imin = lbound(bx);
      const Dim3 imax = ubound(bx);

      constexpr int di = 1;
      const int dj = di * (amax.x - amin.x + 1);
      const int dk = dj * (amax.y - amin.y + 1);

      const Array4<CCTK_REAL> &vars = ghext->mfab.array(mfi);
      CCTK_REAL *restrict const phi = vars.ptr(0, 0, 0, 0);

      for (int k = imin.z; k <= imax.z; ++k)
        for (int j = imin.y; j <= imax.y; ++j)
          for (int i = imin.x; i <= imax.x; ++i) {
            const int idx = i * di + j * dj + k * dk;
            phi[idx] = 0;
            ++count;
          }
    }
    ParallelDescriptor::ReduceIntSum(count);
    assert(count == ghext->ncells * ghext->ncells * ghext->ncells);
  }
  ghext->mfab.FillBoundary(ghext->geom.periodicity());

#pragma omp parallel
  for (MFIter mfi(ghext->mfab,
                  MFItInfo().SetDynamic(true).EnableTiling({1024000, 16, 32}));
       mfi.isValid(); ++mfi) {
    const Box &fbx = mfi.fabbox();
    const Box &bx = mfi.growntilebox();

    const Dim3 amin = lbound(fbx);
    const Dim3 amax = ubound(fbx);
    const Dim3 imin = lbound(bx);
    const Dim3 imax = ubound(bx);

    constexpr int di = 1;
    const int dj = di * (amax.x - amin.x + 1);
    const int dk = dj * (amax.y - amin.y + 1);

    const Array4<CCTK_REAL> &vars = ghext->mfab.array(mfi);
    atomic<CCTK_REAL> *restrict const phi =
        (atomic<CCTK_REAL> *)vars.ptr(0, 0, 0, 0);

    for (int k = imin.z; k <= imax.z; ++k)
      for (int j = imin.y; j <= imax.y; ++j)
#pragma omp simd
        for (int i = imin.x; i <= imax.x; ++i) {
          const int idx = i * di + j * dj + k * dk;
          // phi[idx] += 1.0;
          CCTK_REAL expected = 0.0;
          bool success = phi[idx].compare_exchange_strong(expected, 1.0);
          // assert(success);
        }
  }

  for (MFIter mfi(ghext->mfab); mfi.isValid(); ++mfi) {
    const Box &fbx = mfi.fabbox();
    const Box &bx = mfi.fabbox();

    const Dim3 amin = lbound(fbx);
    const Dim3 amax = ubound(fbx);
    const Dim3 imin = lbound(bx);
    const Dim3 imax = ubound(bx);

    constexpr int di = 1;
    const int dj = di * (amax.x - amin.x + 1);
    const int dk = dj * (amax.y - amin.y + 1);

    const Array4<CCTK_REAL> &vars = ghext->mfab.array(mfi);
    CCTK_REAL *restrict const phi = vars.ptr(0, 0, 0, 0);

    for (int k = imin.z; k <= imax.z; ++k)
      for (int j = imin.y; j <= imax.y; ++j)
        for (int i = imin.x; i <= imax.x; ++i) {
          const int idx = i * di + j * dj + k * dk;
          assert(phi[idx] == 1);
        }
  }
}

extern "C" void WaveToyAMReX_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_VINFO("Initialize iteration %d", cctk_iteration);

  ghext->time = 0.0;
  ghext->delta_time = 0.5 / ghext->ncells;

  const CCTK_REAL t0 = ghext->time;
  const CCTK_REAL dt = ghext->delta_time;
  const CCTK_REAL *restrict const x0 = ghext->geom.ProbLo();
  const CCTK_REAL *restrict const x1 = ghext->geom.ProbHi();

// Initialize phi
#pragma omp parallel
  for (MFIter mfi(ghext->mfab,
                  MFItInfo().SetDynamic(true).EnableTiling({1024000, 16, 32}));
       mfi.isValid(); ++mfi) {
    const Box &fbx = mfi.fabbox();
    const Box &bx = mfi.growntilebox();

    const Dim3 amin = lbound(fbx);
    const Dim3 amax = ubound(fbx);
    const Dim3 imin = lbound(bx);
    const Dim3 imax = ubound(bx);

    constexpr int di = 1;
    const int dj = di * (amax.x - amin.x + 1);
    const int dk = dj * (amax.y - amin.y + 1);

    const Array4<CCTK_REAL> &vars = ghext->mfab.array(mfi);
    CCTK_REAL *restrict const phi = vars.ptr(0, 0, 0, 0);
    CCTK_REAL *restrict const phi_p = vars.ptr(0, 0, 0, 1);

    for (int k = imin.z; k <= imax.z; ++k)
      for (int j = imin.y; j <= imax.y; ++j)
#pragma omp simd
        for (int i = imin.x; i <= imax.x; ++i) {
          const int idx = i * di + j * dj + k * dk;
          CCTK_REAL x = linterp(x0[0], x1[0], -1, 2 * ghext->ncells - 1, 2 * i);
          CCTK_REAL y = linterp(x0[1], x1[1], -1, 2 * ghext->ncells - 1, 2 * j);
          CCTK_REAL z = linterp(x0[2], x1[2], -1, 2 * ghext->ncells - 1, 2 * k);
          phi[idx] = standing(t0, x, y, z);
          phi_p[idx] = standing(t0 - dt, x, y, z);
        }
  }
}

extern "C" void WaveToyAMReX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_VINFO("Evolve iteration %d", cctk_iteration);

  // Cycle time levels
  MultiFab::Copy(ghext->mfab, ghext->mfab, 1, 2, 1, ghext->nghostzones);
  MultiFab::Copy(ghext->mfab, ghext->mfab, 0, 1, 1, ghext->nghostzones);

  const CCTK_REAL dt = ghext->delta_time;
  const CCTK_REAL *restrict const dx = ghext->geom.CellSize();

  // Evolve phi
#pragma omp parallel
  for (MFIter mfi(ghext->mfab,
                  MFItInfo().SetDynamic(true).EnableTiling({1024000, 16, 32}));
       mfi.isValid(); ++mfi) {
    const Box &fbx = mfi.fabbox();
    const Box &bx = mfi.tilebox();

    const Dim3 amin = lbound(fbx);
    const Dim3 amax = ubound(fbx);
    const Dim3 imin = lbound(bx);
    const Dim3 imax = ubound(bx);

    constexpr int di = 1;
    const int dj = di * (amax.x - amin.x + 1);
    const int dk = dj * (amax.y - amin.y + 1);

    const Array4<CCTK_REAL> &vars = ghext->mfab.array(mfi);
    CCTK_REAL *restrict const phi = vars.ptr(0, 0, 0, 0);
    const CCTK_REAL *restrict const phi_p = vars.ptr(0, 0, 0, 1);
    const CCTK_REAL *restrict const phi_p_p = vars.ptr(0, 0, 0, 2);

    for (int k = imin.z; k <= imax.z; ++k)
      for (int j = imin.y; j <= imax.y; ++j)
#pragma omp simd
        for (int i = imin.x; i <= imax.x; ++i) {
          const int idx = i * di + j * dj + k * dk;
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

  // Synchronize
  ghext->mfab.FillBoundary(ghext->geom.periodicity());
}

extern "C" void WaveToyAMReX_Error(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_VINFO("Error iteration %d", cctk_iteration);

  const CCTK_REAL t0 = ghext->time;
  const CCTK_REAL *restrict const x0 = ghext->geom.ProbLo();
  const CCTK_REAL *restrict const x1 = ghext->geom.ProbHi();

// Initialize phi
#pragma omp parallel
  for (MFIter mfi(ghext->mfab,
                  MFItInfo().SetDynamic(true).EnableTiling({1024000, 16, 32}));
       mfi.isValid(); ++mfi) {
    const Box &fbx = mfi.fabbox();
    const Box &bx = mfi.growntilebox();

    const Dim3 amin = lbound(fbx);
    const Dim3 amax = ubound(fbx);
    const Dim3 imin = lbound(bx);
    const Dim3 imax = ubound(bx);

    constexpr int di = 1;
    const int dj = di * (amax.x - amin.x + 1);
    const int dk = dj * (amax.y - amin.y + 1);

    const Array4<CCTK_REAL> &vars = ghext->mfab.array(mfi);
    CCTK_REAL *restrict const err = vars.ptr(0, 0, 0, 3);
    const CCTK_REAL *restrict const phi = vars.ptr(0, 0, 0, 0);

    for (int k = imin.z; k <= imax.z; ++k)
      for (int j = imin.y; j <= imax.y; ++j)
#pragma omp simd
        for (int i = imin.x; i <= imax.x; ++i) {
          const int idx = i * di + j * dj + k * dk;
          CCTK_REAL x = linterp(x0[0], x1[0], -1, 2 * ghext->ncells - 1, 2 * i);
          CCTK_REAL y = linterp(x0[1], x1[1], -1, 2 * ghext->ncells - 1, 2 * j);
          CCTK_REAL z = linterp(x0[2], x1[2], -1, 2 * ghext->ncells - 1, 2 * k);
          err[idx] = phi[idx] - standing(t0, x, y, z);
        }
  }
}

extern "C" void WaveToyAMReX_Output(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_VINFO("Output iteration %d", cctk_iteration);

  // Output phi
  string filename = amrex::Concatenate("wavetoy/phi", cctk_iteration, 6);
  WriteSingleLevelPlotfile(filename, ghext->mfab,
                           {"phi", "phi_p", "phi_p_p", "error"}, ghext->geom,
                           ghext->time, cctk_iteration);
}

} // namespace WaveToyAMReX
