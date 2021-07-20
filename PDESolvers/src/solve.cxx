#include "pdesolvers.hxx"

// TODO: Don't include files from other thorn; create a proper interface
#include "../../CarpetX/src/driver.hxx"
#include "../../CarpetX/src/reduction.hxx"
#include "../../CarpetX/src/schedule.hxx"

#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <AMReX_MultiFabUtil.H>

#include <mpi.h>

#include <petscdm.h>
#include <petscsnes.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <optional>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace PDESolvers {

const int tl = 0;

// TODO: Generalize this
const Arith::vect<int, 3> indextype{0, 0, 0}; // vertex centred

// Tuning parameters for matrix preallocation
// const int dnz = 7;
// const int onz = 4;
const int dnz = 13;
const int onz = 8;

////////////////////////////////////////////////////////////////////////////////

// Erik Schnetter the blocks for vertex centred grids overlap by one
// point at the block boundaries (see
// https://amrex-codes.github.io/amrex/docs_html/Basics.html#id4).
// does FillPatchSingleLevel ensure these points are consistent? or is
// the domain double-valued there, and only the ghost points (not
// shown in this figure) are set?

// weiqun If the original source data are consistent, the results of
// FillPatchSingleLevel are also consistent.  Otherwise, it's undefined which
// original values will used for which boxes.  It depends on how boxes are
// distributed.

// If you want to the data consistent, you can use MultiFab's OverrideSync
// functions.
//     void OverrideSync (const Periodicity& period =
//         Periodicity::NonPeriodic());
//     void OverrideSync (const iMultiFab& msk,
//         const Periodicity& period = Periodicity::NonPeriodic());

// The firs version will use data from the grid with the lowest grid number to
// override others.

// In the second version, you can pass your own mask.

// Erik Schnetter thanks!

// weiqun Even if you want to use the first version's mask, you might still want
// to call
//     //! Owner is the grid with the lowest grid number containing the data.
//     std::unique_ptr<iMultiFab> OwnerMask (const Periodicity& period =
//         Periodicity::NonPeriodic()) const;
// and then use the second function.

// If you need to call OverrideSync repeatedly, you could save the mask.

// Erik Schnetter thanks!

////////////////////////////////////////////////////////////////////////////////

// 1: interior (defined by user)
// 2: outer boundary (also defined by users)
// 3: defined via synchronization (copied from another box on same level)
// 4: defined via prolongation (interpolated from next coarser level)
// 5: defined via restriction (injected from next finer level)
enum class point_type_t {
  undf = 0,
  intr = 1,
  bdry = 2,
  sync = 3,
  prol = 4,
  rest = 5
};

void define_point_type(const cGH *const cctkGH) {
  // Decode Cactus variables
  const int vn_pt = CCTK_VarIndex("PDESolvers::point_type");
  assert(vn_pt >= 0);
  const int gi_pt = CCTK_GroupIndexFromVarI(vn_pt);
  assert(gi_pt >= 0);
  const int v0_pt = CCTK_FirstVarIndexI(gi_pt);
  assert(v0_pt >= 0);
  const int vi_pt = vn_pt - v0_pt;
  assert(vi_pt >= 0);
  const int vn_ind = CCTK_VarIndex("PDESolvers::communication_indicator");
  assert(vn_ind >= 0);
  const int gi_ind = CCTK_GroupIndexFromVarI(vn_ind);
  assert(gi_ind >= 0);
  const int v0_ind = CCTK_FirstVarIndexI(gi_ind);
  assert(v0_ind >= 0);
  const int vi_ind = vn_ind - v0_ind;
  assert(vi_ind >= 0);

  // Initialize point type everywhere, assuming there is no synchronization,
  // prolongation, or restriction
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_pt(
        layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             gf_pt(p.I) = int(point_type_t::undf);
                           });
    grid.loop_int<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             assert(gf_pt(p.I) == int(point_type_t::undf));
                             gf_pt(p.I) = int(point_type_t::intr);
                           });
    grid.loop_bnd<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             assert(gf_pt(p.I) == int(point_type_t::undf));
                             gf_pt(p.I) = int(point_type_t::bdry);
                           });
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    leveldata.groupdata.at(gi_pt)->valid.at(tl).at(vi_pt).set(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::define_point_type"; });
  });

  // Set indicator to level
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_ind(
        layout1,
        static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p)
                               ARITH_INLINE { gf_ind(p.I) = level; });
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    leveldata.groupdata.at(gi_ind)->valid.at(tl).at(vi_ind).set(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::define_point_type before restricting"; });
  });

  // Restrict
  for (int level = int(CarpetX::ghext->leveldata.size()) - 2; level >= 0;
       --level) {
    auto &leveldata = CarpetX::ghext->leveldata.at(level);
    const auto &fineleveldata = CarpetX::ghext->leveldata.at(level + 1);
    amrex::MultiFab &mfab_ind = *leveldata.groupdata.at(gi_ind)->mfab.at(tl);
    const amrex::MultiFab &finemfab_ind =
        *fineleveldata.groupdata.at(gi_ind)->mfab.at(tl);
    const amrex::IntVect reffact{2, 2, 2};
    const int rank = sum(indextype);
    switch (rank) {
    case 0:
      average_down_nodal(finemfab_ind, mfab_ind, reffact);
      break;
    case 1:
      average_down_edges(finemfab_ind, mfab_ind, reffact);
      break;
    case 2:
      average_down_faces(finemfab_ind, mfab_ind, reffact);
      break;
    case 3:
      average_down(finemfab_ind, mfab_ind, 0, 1 /*nvars*/, reffact);
      break;
    default:
      assert(0);
    }
  }

  // Check where the indicator changed; these are the restricted points
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_pt(
        layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const Loop::GF3D2<const CCTK_REAL> gf_ind(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             if (gf_ind(p.I) != level) {
                               assert(false); // TODO: just for testing
                               assert(gf_pt(p.I) != int(point_type_t::bdry));
                               gf_pt(p.I) = int(point_type_t::rest);
                             }
                           });
  });

  // Set indicator to index
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_ind(
        layout1,
        static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p)
                               ARITH_INLINE { gf_ind(p.I) = index; });
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    leveldata.groupdata.at(gi_ind)->valid.at(tl).at(vi_ind).set(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::define_point_type before synchronizing"; });
  });

  // Synchronize (as if there was no prolongation)
  for (int level = 0; level < int(CarpetX::ghext->leveldata.size()); ++level) {
    auto &leveldata = CarpetX::ghext->leveldata.at(level);
    amrex::MultiFab &mfab_ind = *leveldata.groupdata.at(gi_ind)->mfab.at(tl);
    auto physbc_bcs = get_boundaries(*leveldata.groupdata.at(gi_ind));
    CarpetX::CarpetXPhysBCFunct &physbc = std::get<0>(physbc_bcs);
    FillPatchSingleLevel(mfab_ind, 0.0, {&mfab_ind}, {0.0}, 0, 0, 1 /*nvars*/,
                         CarpetX::ghext->amrcore->Geom(level), physbc, 0);
  }

  // Check where the indicator changed; these are the synchronized points.
  // If points are both restricted and synchronized, then we count them as
  // synchronized.
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_pt(
        layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const Loop::GF3D2<const CCTK_REAL> gf_ind(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             if (gf_ind(p.I) != index) {
                               assert(gf_pt(p.I) == int(point_type_t::undf) ||
                                      gf_pt(p.I) == int(point_type_t::bdry));
                               gf_pt(p.I) = int(point_type_t::sync);
                             }
                             // Because we don't support
                             // prolongation/restriction yet
                             if (gf_ind(p.I) == index)
                               assert(gf_pt(p.I) == int(point_type_t::intr) ||
                                      gf_pt(p.I) == int(point_type_t::bdry));
                           });
  });

  // Set indicator to level
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_ind(
        layout1,
        static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p)
                               ARITH_INLINE { gf_ind(p.I) = level; });
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    leveldata.groupdata.at(gi_ind)->valid.at(tl).at(vi_ind).set(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::define_point_type before prolongating"; });
  });

  // Prolongate and synchronize (we cannot just prolongate)
  for (int level = 1; level < int(CarpetX::ghext->leveldata.size()); ++level) {
    auto &leveldata = CarpetX::ghext->leveldata.at(level);
    const auto &coarseleveldata = CarpetX::ghext->leveldata.at(level - 1);
    amrex::MultiFab &mfab_ind = *leveldata.groupdata.at(gi_ind)->mfab.at(tl);
    amrex::MultiFab &coarsemfab_ind =
        *coarseleveldata.groupdata.at(gi_ind)->mfab.at(tl);
    const amrex::IntVect reffact{2, 2, 2};
    amrex::Interpolater *const interpolator =
        CarpetX::get_interpolator(std::array<int, 3>(indextype));
    auto physbc_bcs = get_boundaries(*leveldata.groupdata.at(gi_ind));
    const amrex::Vector<amrex::BCRec> &bcs = std::get<1>(physbc_bcs);
    CarpetX::CarpetXPhysBCFunct &physbc = std::get<0>(physbc_bcs);
    FillPatchTwoLevels(mfab_ind, 0.0, {&coarsemfab_ind}, {0.0}, {&mfab_ind},
                       {0.0}, 0, 0, 1 /*nvars*/,
                       CarpetX::ghext->amrcore->Geom(level - 1),
                       CarpetX::ghext->amrcore->Geom(level), physbc, 0, physbc,
                       0, reffact, interpolator, bcs, 0);
  }

  // Check where the indicator changed; these are the prolongated points. If
  // points are both restricted and prolongated, then we count them as
  // prolongated.
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_pt(
        layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const Loop::GF3D2<const CCTK_REAL> gf_ind(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_ind)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
          if (gf_ind(p.I) != level && gf_pt(p.I) != int(point_type_t::sync)) {
            assert(false); // TODO: just for testing
            assert(gf_pt(p.I) != int(point_type_t::bdry));
            gf_pt(p.I) = int(point_type_t::prol);
          }
        });
  });

  // Invalidate indicator
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    leveldata.groupdata.at(gi_ind)->valid.at(tl).at(vi_ind).set(
        CarpetX::valid_t(), []() { return "PDESolver::define_point_type"; });
  });

  // Collect some statistics
  Arith::vect<int, 6> npoints{0, 0, 0, 0, 0, 0};
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const CarpetX::GridDescBase grid(cctkGH);
    Arith::vect<int, 6> npoints1{0, 0, 0, 0, 0, 0};
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             const int pt = int(gf_pt(p.I));
                             assert(pt >= 0 && pt < 6);
                             ++npoints1[pt];
                           });
    for (int n = 0; n < 6; ++n)
#pragma omp atomic
      npoints[n] += npoints1[n];
  });
  CCTK_VINFO("Point type counts:");
  CCTK_VINFO("  undefined:    %d", npoints[int(point_type_t::undf)]);
  CCTK_VINFO("  interior:     %d", npoints[int(point_type_t::intr)]);
  CCTK_VINFO("  boundary:     %d", npoints[int(point_type_t::bdry)]);
  CCTK_VINFO("  synchronized: %d", npoints[int(point_type_t::sync)]);
  CCTK_VINFO("  prolongated:  %d", npoints[int(point_type_t::prol)]);
  CCTK_VINFO("  restricted:   %d", npoints[int(point_type_t::rest)]);
  CCTK_VINFO("  total:        %d", sum(npoints));
  assert(npoints[int(point_type_t::undf)] == 0);
}

////////////////////////////////////////////////////////////////////////////////

void enumerate_points(const cGH *const cctkGH, int &restrict npoints_local,
                      int &restrict npoints_global,
                      std::vector<std::vector<int> > &restrict block_offsets,
                      std::vector<std::vector<int> > &restrict block_sizes) {
  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  static const int vn_idx = CCTK_VarIndex("PDESolvers::idx");
  assert(vn_idx >= 0);
  static const int gi_idx = CCTK_GroupIndexFromVarI(vn_idx);
  assert(gi_idx >= 0);
  static const int v0_idx = CCTK_FirstVarIndexI(gi_idx);
  assert(v0_idx >= 0);
  static const int vi_idx = vn_idx - v0_idx;
  assert(vi_idx >= 0);
  const int vn_pt = CCTK_VarIndex("PDESolvers::point_type");
  assert(vn_pt >= 0);
  const int gi_pt = CCTK_GroupIndexFromVarI(vn_pt);
  assert(gi_pt >= 0);
  const int v0_pt = CCTK_FirstVarIndexI(gi_pt);
  assert(v0_pt >= 0);
  const int vi_pt = vn_pt - v0_pt;
  assert(vi_pt >= 0);

  // Determine number of blocks per level
  std::vector<int> level_sizes(CarpetX::ghext->leveldata.size(), 0);
  std::vector<int> level_maxblocks(CarpetX::ghext->leveldata.size(), -1);
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
#pragma omp atomic
    ++level_sizes.at(level);
    int &maxblock = level_maxblocks.at(level);
    using std::max;
#pragma omp critical
    maxblock = max(maxblock, block);
  });

  // Allocate data structure
  block_offsets.resize(CarpetX::ghext->leveldata.size()); // process local
  block_sizes.resize(CarpetX::ghext->leveldata.size());
  for (std::size_t level = 0; level < CarpetX::ghext->leveldata.size();
       ++level) {
    assert(level_maxblocks.at(level) + 1 == level_sizes.at(level));
    block_offsets.at(level).resize(level_sizes.at(level));
    block_sizes.at(level).resize(level_sizes.at(level));
  }

  // Enumerate and count points
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const CarpetX::GridDescBase grid(cctkGH);
    int npoints = 0;
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             switch (int(gf_pt(p.I))) {
                             case int(point_type_t::intr):
                             case int(point_type_t::bdry):
                             case int(point_type_t::rest):
                               ++npoints;
                               break;
                             case int(point_type_t::sync):
                             case int(point_type_t::prol):
                               // ignore this point
                               break;
                             default:
                               assert(0);
                             }
                           });
    block_sizes.at(level).at(block) = npoints;
  });

  // Local exclusive prefix sum
  int npoints = 0;
  for (std::size_t level = 0; level < block_offsets.size(); ++level) {
    for (std::size_t block = 0; block < block_offsets.at(level).size();
         ++block) {
      block_offsets.at(level).at(block) = npoints;
      npoints += block_sizes.at(level).at(block);
    }
  }
  npoints_local = npoints;

  // Global exclusive prefix sum
  int npoints_offset = 0;
  MPI_Exscan(&npoints_local, &npoints_offset, 1, MPI_INT, MPI_SUM,
             MPI_COMM_WORLD);
  MPI_Allreduce(&npoints_local, &npoints_global, 1, MPI_INT, MPI_SUM,
                MPI_COMM_WORLD);

  // Set indices
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<CCTK_REAL> gf_idx(
        layout1,
        static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_idx)));
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const CarpetX::GridDescBase grid(cctkGH);
    int idx = npoints_offset + block_offsets.at(level).at(block);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             switch (int(gf_pt(p.I))) {
                             case int(point_type_t::intr):
                             case int(point_type_t::bdry):
                             case int(point_type_t::rest):
                               // solve for this point
                               gf_idx(p.I) = idx++;
                               break;
                             case int(point_type_t::sync):
                             case int(point_type_t::prol):
                               // ignore this point
                               gf_idx(p.I) = -1;
                               break;
                             default:
                               assert(0);
                             }
                           });
    assert(idx == npoints_offset + block_offsets.at(level).at(block) +
                      block_sizes.at(level).at(block));
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    leveldata.groupdata.at(gi_idx)->valid.at(tl).at(vi_idx).set(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::enumerate_points"; });
  });

  // Synchronize index
  for (int level = 0; level < int(CarpetX::ghext->leveldata.size()); ++level) {
    auto &leveldata = CarpetX::ghext->leveldata.at(level);
    amrex::MultiFab &mfab_idx = *leveldata.groupdata.at(gi_idx)->mfab.at(tl);
    auto physbc_bcs = get_boundaries(*leveldata.groupdata.at(gi_idx));
    CarpetX::CarpetXPhysBCFunct &physbc = std::get<0>(physbc_bcs);
    FillPatchSingleLevel(mfab_idx, 0.0, {&mfab_idx}, {0.0}, 0, 0, 1 /*nvars*/,
                         CarpetX::ghext->amrcore->Geom(level), physbc, 0);
  }

  // Check that synchronization worked
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    const Loop::GF3D2<const CCTK_REAL> gf_idx(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_idx)));
    const CarpetX::GridDescBase grid(cctkGH);
    grid.loop_all<0, 0, 0>(grid.nghostzones,
                           [&](const Loop::PointDesc &p) ARITH_INLINE {
                             switch (int(gf_pt(p.I))) {
                             case int(point_type_t::intr):
                             case int(point_type_t::bdry):
                             case int(point_type_t::sync):
                             case int(point_type_t::rest):
                               assert(gf_idx(p.I) >= 0);
                               break;
                             case int(point_type_t::prol):
                               assert(false); // not yet supported
                               break;
                             default:
                               assert(0);
                             }
                           });
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    leveldata.groupdata.at(gi_idx)->valid.at(tl).at(vi_idx).set(
        CarpetX::make_valid_all(),
        []() { return "PDESolver::enumerate_points"; });
  });
}

void copy_Cactus_to_PETSc(const cGH *const cctkGH, Vec vec,
                          const std::vector<int> &varinds,
                          const std::vector<std::vector<int> > &block_offsets,
                          const std::vector<std::vector<int> > &block_sizes) {
  PetscErrorCode ierr;

  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  const int nvars = varinds.size();
  std::vector<int> gis, vis;
  gis.reserve(nvars);
  vis.reserve(nvars);
  for (int vn : varinds) {
    assert(vn >= 0);
    const int gi = CCTK_GroupIndexFromVarI(vn);
    assert(gi >= 0);
    gis.push_back(gi);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = vn - v0;
    assert(vi >= 0);
    vis.push_back(vi);
  }
  static const int vn_pt = CCTK_VarIndex("PDESolvers::point_type");
  assert(vn_pt >= 0);
  static const int gi_pt = CCTK_GroupIndexFromVarI(vn_pt);
  assert(gi_pt >= 0);
  static const int v0_pt = CCTK_FirstVarIndexI(gi_pt);
  assert(v0_pt >= 0);
  static const int vi_pt = vn_pt - v0_pt;
  assert(vi_pt >= 0);

  // Get vector from PETSc
  CCTK_REAL *vec_ptr;
  ierr = VecGetArray(vec, &vec_ptr);
  assert(!ierr);
  CCTK_REAL *restrict const petsc_ptr = vec_ptr;

  // Copy Cactus vector to PETSc
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const CarpetX::GridDescBase grid(cctkGH);
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    for (int n = 0; n < nvars; ++n)
      CarpetX::error_if_invalid(*leveldata.groupdata.at(gis.at(n)), vis.at(n),
                                tl, CarpetX::make_valid_all(), []() {
                                  return "PDESolver::copy_Cactus_to_PETSc";
                                });
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    std::vector<Loop::GF3D2<const CCTK_REAL> > gfs;
    gfs.reserve(nvars);
    for (int n = 0; n < nvars; ++n)
      gfs.emplace_back(layout1, static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(
                                    cctkGH, tl, varinds.at(n))));
    const int block_offset = block_offsets.at(level).at(block);
    int nelems = 0;
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
          switch (int(gf_pt(p.I))) {
          case int(point_type_t::intr):
          case int(point_type_t::bdry):
            for (int n = 0; n < nvars; ++n)
              petsc_ptr[nvars * block_offset + nelems++] = gfs.at(n)(p.I);
            break;
          case int(point_type_t::rest):
          case int(point_type_t::sync):
          case int(point_type_t::prol):
            // ignore this point
            break;
          default:
            assert(0);
          }
        });
    assert(nelems == nvars * block_sizes.at(level).at(block));
  });

  ierr = VecRestoreArray(vec, &vec_ptr);
  assert(!ierr);
}

void copy_PETSc_to_Cactus(const cGH *const cctkGH, Vec vec,
                          const std::vector<int> &varinds,
                          const std::vector<std::vector<int> > &block_offsets,
                          const std::vector<std::vector<int> > &block_sizes) {
  PetscErrorCode ierr;

  int myproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &myproc);

  // Decode Cactus variables
  const int nvars = varinds.size();
  std::vector<int> gis, vis;
  gis.reserve(nvars);
  vis.reserve(nvars);
  for (int vn : varinds) {
    assert(vn >= 0);
    const int gi = CCTK_GroupIndexFromVarI(vn);
    assert(gi >= 0);
    gis.push_back(gi);
    const int v0 = CCTK_FirstVarIndexI(gi);
    assert(v0 >= 0);
    const int vi = vn - v0;
    assert(vi >= 0);
    vis.push_back(vi);
  }
  static const int vn_pt = CCTK_VarIndex("PDESolvers::point_type");
  assert(vn_pt >= 0);
  static const int gi_pt = CCTK_GroupIndexFromVarI(vn_pt);
  assert(gi_pt >= 0);
  static const int v0_pt = CCTK_FirstVarIndexI(gi_pt);
  assert(v0_pt >= 0);
  static const int vi_pt = vn_pt - v0_pt;
  assert(vi_pt >= 0);

  // Get vector from PETSc
  const CCTK_REAL *vec_ptr;
  ierr = VecGetArrayRead(vec, &vec_ptr);
  assert(!ierr);
  const CCTK_REAL *restrict const petsc_ptr = vec_ptr;

  // Copy PETSc vector to Cactus
  CarpetX::loop_over_blocks(cctkGH, [&](const int level, const int index,
                                        const int block,
                                        const cGH *restrict const cctkGH) {
    const Loop::GF3D2layout layout1(cctkGH, indextype);
    const CarpetX::GridDescBase grid(cctkGH);
    const auto &leveldata = CarpetX::ghext->leveldata.at(level);
    const Loop::GF3D2<const CCTK_REAL> gf_pt(
        layout1,
        static_cast<const CCTK_REAL *>(CCTK_VarDataPtrI(cctkGH, tl, vn_pt)));
    std::vector<Loop::GF3D2<CCTK_REAL> > gfs;
    gfs.reserve(nvars);
    for (int n = 0; n < nvars; ++n)
      gfs.emplace_back(layout1, static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                                    cctkGH, tl, varinds.at(n))));
    const int block_offset = block_offsets.at(level).at(block);
    int nelems = 0;
    grid.loop_all<0, 0, 0>(
        grid.nghostzones, [&](const Loop::PointDesc &p) ARITH_INLINE {
          switch (int(gf_pt(p.I))) {
          case int(point_type_t::intr):
          case int(point_type_t::bdry):
            for (int n = 0; n < nvars; ++n)
              gfs.at(n)(p.I) = petsc_ptr[nvars * block_offset + nelems++];
            break;
          case int(point_type_t::rest):
          case int(point_type_t::sync):
          case int(point_type_t::prol):
            // ignore this point
            break;
          default:
            assert(0);
          }
        });
    assert(nelems == nvars * block_sizes.at(level).at(block));
    for (int n = 0; n < nvars; ++n)
      leveldata.groupdata.at(gis.at(n))->valid.at(tl).at(vis.at(n)).set(
          CarpetX::make_valid_int() | CarpetX::make_valid_outer(),
          []() { return "PDESolver::copy_PETSc_to_Cactus"; });
  });

  ierr = VecRestoreArrayRead(vec, &vec_ptr);
  assert(!ierr);
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void PDESolvers_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Setup;
  DECLARE_CCTK_PARAMETERS;

  // Parse Cactus parameter into PETSc command line options
  /*static*/ std::vector<std::string> args;
  args.push_back("cactus");
  std::istringstream buf(petsc_options);
  std::copy(std::istream_iterator<std::string>(buf),
            std::istream_iterator<std::string>(), std::back_inserter(args));

  /*static*/ std::vector<char *> argptrs;
  for (auto &arg : args)
    argptrs.push_back(arg.data());
  argptrs.push_back(nullptr);

  /*static*/ int argc = args.size();
  /*static*/ char **argv = argptrs.data();
  PetscErrorCode ierr = PetscInitialize(&argc, &argv, NULL, NULL);
  assert(!ierr);
}

namespace {
PetscErrorCode FormFunction(SNES snes, Vec x, Vec f, void *fun) {
  return (
      *static_cast<std::function<PetscErrorCode(SNES snes, Vec x, Vec f)> *>(
          fun))(snes, x, f);
}

PetscErrorCode FormJacobian(SNES snes, Vec x, Mat J, Mat B, void *fun) {
  return (*static_cast<
          std::function<PetscErrorCode(SNES snes, Vec x, Mat J, Mat B)> *>(
      fun))(snes, x, J, B);
}
} // namespace

extern "C" void PDESolvers_Solve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Solve;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("PDESolvers_Solve: Solving on %d levels",
             int(CarpetX::ghext->leveldata.size()));

  PetscErrorCode ierr;

  // Grid layout

  define_point_type(cctkGH);

  int npoints_local, npoints_global;
  // block_offsets are process local
  std::vector<std::vector<int> > block_offsets, block_sizes;
  enumerate_points(cctkGH, npoints_local, npoints_global, block_offsets,
                   block_sizes);

  // TODO: fix this
  const std::vector<int> solinds{CCTK_VarIndex("Poisson2::sol")};
  const std::vector<int> resinds{CCTK_VarIndex("Poisson2::res")};
  const int nvars = solinds.size();
  assert(int(resinds.size()) == nvars);

  // Nonlinear solver

  SNES snes;
  ierr = SNESCreate(PETSC_COMM_WORLD, &snes);
  assert(!ierr);
  ierr = SNESSetFromOptions(snes);
  assert(!ierr);
  ierr = SNESSetAlwaysComputesFinalResidual(snes, PETSC_TRUE);
  assert(!ierr);

  // State vector and evaluation function

  Vec r;
  ierr = VecCreateMPI(PETSC_COMM_WORLD, nvars * npoints_local,
                      nvars * npoints_global, &r);
  assert(!ierr);
  ierr = VecSetFromOptions(r);
  assert(!ierr);

  std::function<PetscErrorCode(SNES snes, Vec x, Vec f)> evalf =
      [&](SNES snes, Vec x, Vec f) {
        copy_PETSc_to_Cactus(cctkGH, x, solinds, block_offsets, block_sizes);
        CallScheduleGroup(cctkGH, "PDESolvers_Residual");
        copy_Cactus_to_PETSc(cctkGH, f, resinds, block_offsets, block_sizes);
        return 0;
      };
  ierr = SNESSetFunction(snes, r, FormFunction, &evalf);
  assert(!ierr);

  // Matrix and Jacobian evaluation function

  Mat J;
  ierr = MatCreateAIJ(PETSC_COMM_WORLD, nvars * npoints_local,
                      nvars * npoints_local, nvars * npoints_global,
                      nvars * npoints_global, dnz, NULL, onz, NULL, &J);
  assert(!ierr);
  ierr = MatSetFromOptions(J);
  assert(!ierr);
  jacobians = std::make_optional<jacobians_t>();

  std::function<PetscErrorCode(SNES snes, Vec x, Mat J, Mat B)> evalJ =
      [&](SNES snes, Vec x, Mat J, Mat B) {
        copy_PETSc_to_Cactus(cctkGH, x, solinds, block_offsets, block_sizes);
        CallScheduleGroup(cctkGH, "PDESolvers_Jacobian");
        jacobians->define_matrix(J);
        jacobians->clear();
        if (false) {
          MatInfo info;
          PetscErrorCode ierr;
          ierr = MatGetInfo(J, MAT_GLOBAL_SUM, &info);
          assert(!ierr);
          CCTK_VINFO("Jacobian info: nz_allocated=%g nz_used=%g nz_unneeded=%g "
                     "memory=%g",
                     double(info.nz_allocated), double(info.nz_used),
                     double(info.nz_unneeded), double(info.memory));
        }
        return 0;
      };
  ierr = SNESSetJacobian(snes, J, J, FormJacobian, &evalJ);
  assert(!ierr);

  {
    double atol, rtol, stol;
    int maxit, maxf;
    ierr = SNESGetTolerances(snes, &atol, &rtol, &stol, &maxit, &maxf);
    assert(!ierr);
    CCTK_VINFO("Settings: atol=%g, rtol=%g, stol=%g, maxit=%d, maxf=%d", atol,
               rtol, stol, maxit, maxf);
  }

  // Initialize solution

  Vec x;
  ierr = VecDuplicate(r, &x);
  assert(!ierr);
  copy_Cactus_to_PETSc(cctkGH, x, solinds, block_offsets, block_sizes);

  // Solve

  CCTK_INFO("Solve");
  ierr = SNESSolve(snes, NULL, x);
  assert(!ierr);

  {
    SNESConvergedReason reason;
    ierr = SNESGetConvergedReason(snes, &reason);
    assert(!ierr);
    switch (reason) {
    case SNES_CONVERGED_FNORM_ABS:
    case SNES_CONVERGED_FNORM_RELATIVE:
    case SNES_CONVERGED_SNORM_RELATIVE:
    case SNES_CONVERGED_ITS:
      // SNES converged; do nothing
      CCTK_INFO("SNES iterations converged");
      break;
    case SNES_DIVERGED_FUNCTION_DOMAIN:
    case SNES_DIVERGED_FUNCTION_COUNT:
    case SNES_DIVERGED_LINEAR_SOLVE:
    case SNES_DIVERGED_FNORM_NAN:
    case SNES_DIVERGED_MAX_IT:
    case SNES_DIVERGED_LINE_SEARCH:
    case SNES_DIVERGED_INNER:
    case SNES_DIVERGED_LOCAL_MIN:
    case SNES_DIVERGED_DTOL:
    case SNES_DIVERGED_JACOBIAN_DOMAIN:
    case SNES_DIVERGED_TR_DELTA:
    case SNES_CONVERGED_ITERATING:
      CCTK_WARN(CCTK_WARN_ALERT, "****************************************");
      CCTK_WARN(CCTK_WARN_ALERT, "SNES iterations diverged");
      CCTK_WARN(CCTK_WARN_ALERT, "****************************************");
      break;
    default:
      CCTK_ERROR("Unknown SNES converged reason");
    }
  }

  {
    int niters;
    ierr = SNESGetIterationNumber(snes, &niters);
    assert(!ierr);
    int liters;
    ierr = SNESGetLinearSolveIterations(snes, &liters);
    assert(!ierr);
    CCTK_VINFO("Iterations: %d Newton, %d linear", niters, liters);
  }

  // Extract solution

  copy_PETSc_to_Cactus(cctkGH, x, solinds, block_offsets, block_sizes);
  CallScheduleGroup(cctkGH, "PDESolvers_Residual");
  {
    const int nvars = resinds.size();
    std::vector<int> gis, vis;
    gis.reserve(nvars);
    vis.reserve(nvars);
    for (int vn : resinds) {
      assert(vn >= 0);
      const int gi = CCTK_GroupIndexFromVarI(vn);
      assert(gi >= 0);
      gis.push_back(gi);
      const int v0 = CCTK_FirstVarIndexI(gi);
      assert(v0 >= 0);
      const int vi = vn - v0;
      assert(vi >= 0);
      vis.push_back(vi);
    }
    const CarpetX::reduction<CCTK_REAL, 3> res_norm =
        CarpetX::reduce(gis.at(0), vis.at(0), tl);
    std::ostringstream buf;
    buf << res_norm;
    CCTK_VINFO("Residual norm: %s", buf.str().c_str());
  }

  // Clean up

  jacobians.reset();
  VecDestroy(&x);
  VecDestroy(&r);
  MatDestroy(&J);
  SNESDestroy(&snes);
  CCTK_INFO("Done.");
}

extern "C" void PDESolvers_Shutdown(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_PDESolvers_Shutdown;
  PetscFinalize();
}

} // namespace PDESolvers
