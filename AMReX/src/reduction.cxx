#include <driver.hxx>
#include <reduction.hxx>

#include <AMReX_Orientation.H>

namespace AMReX {
using namespace amrex;
using namespace std;

namespace {
// Convert a (direction, face) pair to an AMReX Orientation
Orientation orient(int d, int f) {
  return Orientation(d, Orientation::Side(f));
}
} // namespace

template <typename T>
reduction<T> reduce_array(const T *restrict var, const int *restrict ash,
                          const int *restrict lsh) {
  reduction<T> red;
  constexpr int di = 1;
  const int dj = di * ash[0];
  const int dk = dj * ash[1];
  for (int k = 0; k < lsh[2]; ++k) {
    for (int j = 0; j < lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < lsh[0]; ++i) {
        int idx = di * i + dj * j + dk * k;
        red = reduction<T>(red, reduction<T>(var[idx]));
      }
    }
  }
  return red;
}

reduction<CCTK_REAL> reduce(int gi, int vi, int tl) {
  reduction<CCTK_REAL> red;
#pragma omp parallel reduction(reduction : red)
  {
    for (auto &restrict leveldata : ghext->leveldata) {
      auto &restrict groupdata = leveldata.groupdata.at(gi);
      MultiFab &mfab = *groupdata.mfab.at(tl);
      auto mfitinfo =
          MFItInfo().SetDynamic(true).EnableTiling({1024000, 16, 32});
      for (MFIter mfi(mfab, mfitinfo); mfi.isValid(); ++mfi) {
        const Box &fbx = mfi.fabbox();
        const Box &bx = mfi.tilebox();

        const Dim3 imin = lbound(bx);
        int ash[3], lsh[3];
        for (int d = 0; d < dim; ++d)
          ash[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)] + 1;
        for (int d = 0; d < dim; ++d)
          lsh[d] = bx[orient(d, 1)] - bx[orient(d, 0)] + 1;

        const Array4<CCTK_REAL> &vars = groupdata.mfab.at(tl)->array(mfi);
        const CCTK_REAL *restrict var = vars.ptr(imin.x, imin.y, imin.z, vi);

        red += reduce_array(var, ash, lsh);
      }
    }
  }
  return red;
}
} // namespace AMReX
