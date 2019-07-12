#include <AMReX.hxx>
using namespace AMReX;

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <atomic>
using namespace std;

namespace AMReXTest {

extern "C" void AMReXTest_CheckLoops(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_VINFO("Check loop iterators");

  // Set all grid points to zero: Loop over the interior, then
  // synchronize
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

  // Increase each grid point (both interior and ghost zones) by one
  // using atomic operations. This will catch cases where we either
  // omit points or traverse points twice.
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

  // Check all grid points whether they are set to one
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

} // namespace AMReXTest
