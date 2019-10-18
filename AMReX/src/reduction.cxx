#include "driver.hxx"
#include "reduction.hxx"
#include "schedule.hxx"

#include <cctk_Parameters.h>

#include <AMReX_MultiFabUtil.H>
#include <AMReX_Orientation.H>

#include <cstdlib>

namespace AMReX {
using namespace amrex;
using namespace std;

namespace {
template <typename T>
void mpi_reduce_typed(const void *restrict x0, void *restrict y0,
                      const int *restrict length) {
  const T *restrict x = static_cast<const T *>(x0);
  T *restrict y = static_cast<T *>(y0);
#pragma omp simd
  for (int i = 0; i < *length; ++i)
    y[i] += x[i];
}

void mpi_reduce(void *restrict x, void *restrict y, int *restrict length,
                MPI_Datatype *restrict datatype) {
  if (*datatype == MPI_FLOAT)
    return mpi_reduce_typed<float>(x, y, length);
  if (*datatype == MPI_DOUBLE)
    return mpi_reduce_typed<double>(x, y, length);
  if (*datatype == MPI_LONG_DOUBLE)
    return mpi_reduce_typed<long double>(x, y, length);
  abort();
}
} // namespace

MPI_Op reduction_mpi_op() {
  static MPI_Op op = MPI_OP_NULL;
  if (op == MPI_OP_NULL)
    MPI_Op_create(mpi_reduce, 1, &op);
  return op;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
reduction<T> reduce_array(const Geometry &geom,
                          const GHExt::LevelData::GroupData &groupdata,
                          const GridPtrDesc &grid, const T *restrict ptr,
                          const Array4<const int> *restrict finemask_array4) {
  const CCTK_REAL *restrict dx = geom.CellSize();
  CCTK_REAL dV = 1.0;
  for (int d = 0; d < dim; ++d)
    dV *= dx[d];

  reduction<T> red;
  grid.loop_int(groupdata.indextype, [&](const Loop::PointDesc &p) {
    if (!finemask_array4 || !(*finemask_array4)(grid.cactus_offset.x + p.i,
                                                grid.cactus_offset.y + p.j,
                                                grid.cactus_offset.z + p.k))
      red = reduction<T>(red, reduction<T>(dV, ptr[p.idx]));
  });

  return red;
}

reduction<CCTK_REAL> reduce(int gi, int vi, int tl) {
  DECLARE_CCTK_PARAMETERS;

  cGroup group;
  int ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);
  assert(group.grouptype == CCTK_GF);

  reduction<CCTK_REAL> red;
  for (auto &restrict leveldata : ghext->leveldata) {
    const auto &restrict geom = ghext->amrcore->Geom(leveldata.level);
    const auto &restrict groupdata = leveldata.groupdata.at(gi);
    const MultiFab &mfab = *groupdata.mfab.at(tl);
    unique_ptr<iMultiFab> finemask;

    if (!groupdata.valid.at(tl).at(vi).valid_int) {
      CCTK_VWARN(CCTK_WARN_ALERT,
                 "Variable %s has invalid data in the interior of level %d",
                 CCTK_FullVarName(groupdata.firstvarindex + vi),
                 leveldata.level);
    }

    const int fine_level = leveldata.level + 1;
    if (fine_level < int(ghext->leveldata.size())) {
      const auto &restrict fine_leveldata = ghext->leveldata.at(fine_level);
      const auto &restrict fine_groupdata = fine_leveldata.groupdata.at(gi);
      const MultiFab &fine_mfab = *fine_groupdata.mfab.at(tl);

      const IntVect reffact{2, 2, 2};
      finemask = make_unique<iMultiFab>(
          makeFineMask(mfab, fine_mfab.boxArray(), reffact));
    }

    auto mfitinfo = MFItInfo().SetDynamic(true).EnableTiling(
        {max_tile_size_x, max_tile_size_y, max_tile_size_z});
#pragma omp parallel reduction(reduction : red)
    for (MFIter mfi(*leveldata.mfab0, mfitinfo); mfi.isValid(); ++mfi) {
      GridPtrDesc grid(leveldata, mfi);
      unique_ptr<Array4<const int> > finemask_array4;
      if (finemask)
        finemask_array4 = make_unique<Array4<const int> >(finemask->array(mfi));
      const Array4<const CCTK_REAL> &vars = mfab.array(mfi);
      const CCTK_REAL *restrict ptr = grid.ptr(vars, vi);
      red += reduce_array(geom, groupdata, grid, ptr, finemask_array4.get());
    }
  }

  MPI_Op op = reduction_mpi_op();
  MPI_Allreduce(MPI_IN_PLACE, &red, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);

  return red;
}
} // namespace AMReX
