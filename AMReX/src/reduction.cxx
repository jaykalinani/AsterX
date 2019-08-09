#include "driver.hxx"
#include "reduction.hxx"
#include "schedule.hxx"

#include <cctk_Parameters.h>

#include <AMReX_Orientation.H>

namespace AMReX {
using namespace amrex;
using namespace std;

template <typename T>
reduction<T> reduce_array(const Geometry &geom,
                          const GHExt::LevelData::GroupData &groupdata,
                          const T *restrict ptr, const GridDesc &grid) {
  const CCTK_REAL *restrict dx = geom.CellSize();
  CCTK_REAL dV = 1.0;
  for (int d = 0; d < dim; ++d)
    dV *= dx[d];

  reduction<T> red;
  grid.loop_int(groupdata.indextype, [&](const Loop::PointDesc &p) {
    red = reduction<T>(red, reduction<T>(dV, ptr[p.idx]));
  });

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
      auto mfitinfo = MFItInfo().SetDynamic(true).EnableTiling(
          {max_tile_size_x, max_tile_size_y, max_tile_size_z});
      for (MFIter mfi(*leveldata.mfab0, mfitinfo); mfi.isValid(); ++mfi) {
        GridPtrDesc grid(leveldata, mfi);

        const Array4<CCTK_REAL> &vars = groupdata.mfab.at(tl)->array(mfi);
        const CCTK_REAL *restrict ptr = grid.ptr(vars, vi);

#warning "TODO: Skip refined regions"
        red += reduce_array(geom, groupdata, ptr, grid);
      }
    }
  }
  return red;
}
} // namespace AMReX
