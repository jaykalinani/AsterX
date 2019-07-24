#include <driver.hxx>
#include <reduction.hxx>

#include <cctk_Parameters.h>

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
reduction<T> reduce_array(const Geometry &geom, const T *restrict ptr,
                          const int *restrict ash, const int *restrict lsh) {
  const CCTK_REAL *restrict dx = geom.CellSize();
  CCTK_REAL dV = 1.0;
  for (int d = 0; d < dim; ++d)
    dV *= dx[d];
  reduction<T> red;
  constexpr int di = 1;
  const int dj = di * ash[0];
  const int dk = dj * ash[1];
  for (int k = 0; k < lsh[2]; ++k) {
    for (int j = 0; j < lsh[1]; ++j) {
#pragma omp simd
      for (int i = 0; i < lsh[0]; ++i) {
        int idx = di * i + dj * j + dk * k;
        red = reduction<T>(red, reduction<T>(dV, ptr[idx]));
      }
    }
  }
  return red;
}

reduction<CCTK_REAL> reduce(int gi, int vi, int tl) {
  DECLARE_CCTK_PARAMETERS;

  reduction<CCTK_REAL> red;
#pragma omp parallel reduction(reduction : red)
  {
    for (auto &restrict leveldata : ghext->leveldata) {
      const auto &restrict geom = ghext->amrcore->Geom(leveldata.level);
      auto &restrict groupdata = leveldata.groupdata.at(gi);
      MultiFab &mfab = *groupdata.mfab.at(tl);
      auto mfitinfo = MFItInfo().SetDynamic(true).EnableTiling(
          {max_tile_size_x, max_tile_size_y, max_tile_size_z});
      for (MFIter mfi(mfab, mfitinfo); mfi.isValid(); ++mfi) {
        const Box &fbx = mfi.fabbox();
        const Box &bx = mfi.tilebox();

        const Dim3 imin = lbound(bx);
        int ash[3], lsh[3];
        for (int d = 0; d < dim; ++d)
          ash[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)] + 1;
        // Note: This excludes ghosts, it's not the proper Cactus lsh
        for (int d = 0; d < dim; ++d)
          lsh[d] = bx[orient(d, 1)] - bx[orient(d, 0)] + 1;

        const Array4<CCTK_REAL> &vars = groupdata.mfab.at(tl)->array(mfi);
        const CCTK_REAL *restrict ptr = vars.ptr(imin.x, imin.y, imin.z, vi);

#warning "TODO: Skip refined regions"
        red += reduce_array(geom, ptr, ash, lsh);
      }
    }
  }
  return red;
}
} // namespace AMReX
