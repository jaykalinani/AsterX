#include "driver.hxx"
#include "io.hxx"
#include "prolongate_3d_cc_rf2.hxx"
#include "prolongate_3d_rf2.hxx"
#include "schedule.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <AMReX.H>
#include <AMReX_BCRec.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>

#include <omp.h>
#include <mpi.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <string>

namespace AMReX {
using namespace amrex;
using namespace std;

// Global variables

int ghext_handle = -1;

amrex::AMReX *restrict pamrex = nullptr;
unique_ptr<GHExt> ghext;

// Registered functions

void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cctkGH);
int InitGH(cGH *cctkGH);
int ScheduleTraverseGH(cGH *cctkGH, const char *where);

int MyProc(const cGH *cctkGH);
int nProcs(const cGH *cctkGH);
int Exit(cGH *cctkGH, int retval);
int Abort(cGH *cctkGH, int retval);
int Barrier(const cGH *cctkGHa);

////////////////////////////////////////////////////////////////////////////////

// AmrCore functions

CactusAmrCore::CactusAmrCore() {}
CactusAmrCore::CactusAmrCore(const RealBox *rb, int max_level_in,
                             const Vector<int> &n_cell_in, int coord,
                             Vector<IntVect> ref_ratios, const int *is_per)
    : AmrCore(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per) {}
CactusAmrCore::CactusAmrCore(const RealBox &rb, int max_level_in,
                             const Vector<int> &n_cell_in, int coord,
                             Vector<IntVect> const &ref_ratios,
                             Array<int, AMREX_SPACEDIM> const &is_per)
    : AmrCore(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per) {}

CactusAmrCore::~CactusAmrCore() {}

void CactusAmrCore::ErrorEst(const int level, TagBoxArray &tags, Real time,
                             int ngrow) {
  DECLARE_CCTK_PARAMETERS;

  // Don't regrid before Cactus is ready to
  if (level >= int(ghext->leveldata.size()))
    return;

  if (verbose)
    CCTK_VINFO("ErrorEst level %d", level);

  const int gi = CCTK_GroupIndex("AMReX::regrid_error");
  assert(gi >= 0);
  const int vi = 0;
  const int tl = 0;

  auto &restrict leveldata = ghext->leveldata.at(level);
  auto &restrict groupdata = leveldata.groupdata.at(gi);
  // Ensure the error estimate has been set
  assert(groupdata.valid.at(tl).at(vi).valid_int);
  auto mfitinfo = MFItInfo().SetDynamic(true).EnableTiling(
      {max_tile_size_x, max_tile_size_y, max_tile_size_z});
#pragma omp parallel
  for (MFIter mfi(*leveldata.mfab0, mfitinfo); mfi.isValid(); ++mfi) {
    GridPtrDesc grid(leveldata, mfi);

    const Array4<CCTK_REAL> &vars = groupdata.mfab.at(tl)->array(mfi);
    CCTK_REAL *restrict const err = grid.ptr(vars, vi);
    const Array4<char> &tags_array4 = tags.array(mfi);

    grid.loop_int(groupdata.indextype, [&](const Loop::PointDesc &p) {
      tags_array4(grid.cactus_offset.x + p.i, grid.cactus_offset.y + p.j,
                  grid.cactus_offset.z + p.k) =
          err[p.idx] >= regrid_error_threshold ? TagBox::SET : TagBox::CLEAR;
    });
    grid.loop_bnd(groupdata.indextype, [&](const Loop::PointDesc &p) {
      tags_array4(grid.cactus_offset.x + p.i, grid.cactus_offset.y + p.j,
                  grid.cactus_offset.z + p.k) = TagBox::CLEAR;
    });
  }
}

array<int, dim> get_group_indextype(const int gi) {
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  array<CCTK_INT, dim> index;
  int iret = Util_TableGetIntArray(tags, dim, index.data(), "index");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    index = {1, 1, 1}; // default: cell-centred
  } else if (iret >= 0) {
    assert(iret == dim);
  } else {
    assert(0);
  }
  array<int, dim> indextype;
  for (int d = 0; d < dim; ++d)
    indextype[d] = index[d];
  return indextype;
}

void SetupLevel(int level, const BoxArray &ba, const DistributionMapping &dm) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("SetupLevel level %d", level);

  assert(level == int(ghext->leveldata.size()));
  ghext->leveldata.resize(level + 1);
  GHExt::LevelData &leveldata = ghext->leveldata.at(level);
  leveldata.level = level;
#warning "TODO: Make this an empty MultiFab"
  leveldata.mfab0 = make_unique<MultiFab>(ba, dm, 1, ghost_size);
  assert(ba.ixType() ==
         IndexType(IndexType::CELL, IndexType::CELL, IndexType::CELL));

  const int numgroups = CCTK_NumGroups();
  leveldata.groupdata.resize(numgroups);
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);
    assert(group.grouptype == CCTK_GF);
    assert(group.vartype == CCTK_VARIABLE_REAL);
    assert(group.disttype == CCTK_DISTRIB_DEFAULT);
    assert(group.dim == dim);

    GHExt::LevelData::GroupData &groupdata = leveldata.groupdata.at(gi);
    groupdata.groupindex = gi;
    groupdata.firstvarindex = CCTK_FirstVarIndexI(gi);
    groupdata.numvars = group.numvars;
    groupdata.indextype = get_group_indextype(gi);

    // Allocate grid hierarchies
    const BoxArray &gba = convert(
        ba,
        IndexType(groupdata.indextype[0] ? IndexType::CELL : IndexType::NODE,
                  groupdata.indextype[1] ? IndexType::CELL : IndexType::NODE,
                  groupdata.indextype[2] ? IndexType::CELL : IndexType::NODE));
    groupdata.mfab.resize(group.numtimelevels);
    groupdata.valid.resize(group.numtimelevels);
    for (int tl = 0; tl < int(groupdata.mfab.size()); ++tl) {
      groupdata.mfab.at(tl) =
          make_unique<MultiFab>(gba, dm, groupdata.numvars, ghost_size);
      groupdata.valid.at(tl) = vector<valid_t>(groupdata.numvars);
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        poison_invalid(leveldata, groupdata, vi, tl);
    }
  }
}

Interpolater *get_interpolator() {
  DECLARE_CCTK_PARAMETERS;
  switch (prolongation_order) {
  case 0:
    return &pc_interp;
  case 1:
    return &cell_bilinear_interp;
  case 2:
    return &quadratic_interp;
  case 4:
    return &prolongate3d_cc_rf2_o4;
  }
  assert(0);
}

void CactusAmrCore::MakeNewLevelFromScratch(int level, Real time,
                                            const BoxArray &ba,
                                            const DistributionMapping &dm) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("MakeNewLevelFromScratch level %d", level);

  SetupLevel(level, ba, dm);

  if (saved_cctkGH) {
    assert(current_level == -1);
    current_level = level;
    CCTK_Traverse(saved_cctkGH, "CCTK_BASEGRID");
    // CCTK_Traverse(saved_cctkGH, "CCTK_POSTREGRID");
    current_level = -1;
  }
}

void CactusAmrCore::MakeNewLevelFromCoarse(int level, Real time,
                                           const BoxArray &ba,
                                           const DistributionMapping &dm) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("MakeNewLevelFromCoarse level %d", level);
  assert(level > 0);

  SetupLevel(level, ba, dm);

  // Prolongate
  Interpolater *const interpolator = get_interpolator();
  auto &leveldata = ghext->leveldata.at(level);
  auto &coarseleveldata = ghext->leveldata.at(level - 1);
  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    auto &restrict groupdata = leveldata.groupdata.at(gi);
    auto &restrict coarsegroupdata = coarseleveldata.groupdata.at(gi);
    assert(coarsegroupdata.numvars == groupdata.numvars);
    PhysBCFunctNoOp cphysbc;
    PhysBCFunctNoOp fphysbc;
    const IntVect reffact{2, 2, 2};
    // boundary conditions
    const BCRec bcrec(periodic_x ? BCType::int_dir : BCType::reflect_odd,
                      periodic_y ? BCType::int_dir : BCType::reflect_odd,
                      periodic_z ? BCType::int_dir : BCType::reflect_odd,
                      periodic_x ? BCType::int_dir : BCType::reflect_odd,
                      periodic_y ? BCType::int_dir : BCType::reflect_odd,
                      periodic_z ? BCType::int_dir : BCType::reflect_odd);
    const Vector<BCRec> bcs(groupdata.numvars, bcrec);

    // If there is more than one time level, then we don't prolongate
    // the oldest.
    int ntls = groupdata.mfab.size();
    int prolongate_tl = ntls > 1 ? ntls - 1 : ntls;
    for (int tl = 0; tl < ntls; ++tl)
      groupdata.valid.at(tl) = vector<valid_t>(groupdata.numvars);
    for (int tl = 0; tl < prolongate_tl; ++tl) {
      // Only interpolate if coarse grid data are valid
      bool all_invalid = true;
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        all_invalid &= !coarsegroupdata.valid.at(tl).at(vi).valid_int &&
                       !coarsegroupdata.valid.at(tl).at(vi).valid_bnd;
      if (all_invalid) {
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          groupdata.valid.at(tl).at(vi).valid_int = false;
          groupdata.valid.at(tl).at(vi).valid_bnd = false;
        }
      } else {
        for (int vi = 0; vi < groupdata.numvars; ++vi)
          assert(coarsegroupdata.valid.at(tl).at(vi).valid_int &&
                 coarsegroupdata.valid.at(tl).at(vi).valid_bnd);
        InterpFromCoarseLevel(
            *groupdata.mfab.at(tl), 0.0, *coarsegroupdata.mfab.at(tl), 0, 0,
            groupdata.numvars, ghext->amrcore->Geom(level - 1),
            ghext->amrcore->Geom(level), cphysbc, 0, fphysbc, 0, reffact,
            interpolator, bcs, 0);
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          groupdata.valid.at(tl).at(vi).valid_int =
              coarsegroupdata.valid.at(tl).at(vi).valid_int &&
              coarsegroupdata.valid.at(tl).at(vi).valid_bnd;
          groupdata.valid.at(tl).at(vi).valid_bnd = false;
        }
      }
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        poison_invalid(leveldata, groupdata, vi, tl);
    }
  }

  if (saved_cctkGH) {
    assert(current_level == -1);
    current_level = level;
    CCTK_Traverse(saved_cctkGH, "CCTK_BASEGRID");
    CCTK_Traverse(saved_cctkGH, "CCTK_POSTREGRID");
    current_level = -1;
  }
}

void CactusAmrCore::RemakeLevel(int level, Real time, const BoxArray &ba,
                                const DistributionMapping &dm) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("RemakeLevel level %d", level);

  // Copy or prolongate
  Interpolater *const interpolator = get_interpolator();
  auto &leveldata = ghext->leveldata.at(level);
#warning "TODO: Make this an empty MultiFab"
  leveldata.mfab0 = make_unique<MultiFab>(ba, dm, 1, ghost_size);
  assert(ba.ixType() ==
         IndexType(IndexType::CELL, IndexType::CELL, IndexType::CELL));

  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    auto &restrict groupdata = leveldata.groupdata.at(gi);

    const BoxArray &gba = convert(
        ba,
        IndexType(groupdata.indextype[0] ? IndexType::CELL : IndexType::NODE,
                  groupdata.indextype[1] ? IndexType::CELL : IndexType::NODE,
                  groupdata.indextype[2] ? IndexType::CELL : IndexType::NODE));

    // If there is more than one time level, then we don't
    // prolongate the oldest.
    int ntls = groupdata.mfab.size();
    int prolongate_tl = ntls > 1 ? ntls - 1 : ntls;

    // This must not happen
    assert(leveldata.level > 0);

    // Copy from same level and/or prolongate from next coarser level
    auto &coarseleveldata = ghext->leveldata.at(level - 1);
    auto &restrict coarsegroupdata = coarseleveldata.groupdata.at(gi);
    assert(coarsegroupdata.numvars == groupdata.numvars);

    PhysBCFunctNoOp cphysbc;
    PhysBCFunctNoOp fphysbc;
    const IntVect reffact{2, 2, 2};
    // boundary conditions
    const BCRec bcrec(periodic_x ? BCType::int_dir : BCType::reflect_odd,
                      periodic_y ? BCType::int_dir : BCType::reflect_odd,
                      periodic_z ? BCType::int_dir : BCType::reflect_odd,
                      periodic_x ? BCType::int_dir : BCType::reflect_odd,
                      periodic_y ? BCType::int_dir : BCType::reflect_odd,
                      periodic_z ? BCType::int_dir : BCType::reflect_odd);
    const Vector<BCRec> bcs(groupdata.numvars, bcrec);
    for (int tl = 0; tl < ntls; ++tl) {
      auto mfab = make_unique<MultiFab>(gba, dm, groupdata.numvars, ghost_size);
      auto valid = vector<valid_t>(groupdata.numvars);
      if (poison_undefined_values) {
        // Set new grid functions to nan
        auto mfitinfo = MFItInfo().SetDynamic(true).EnableTiling(
            {max_tile_size_x, max_tile_size_y, max_tile_size_z});
#pragma omp parallel
        for (MFIter mfi(*leveldata.mfab0, mfitinfo); mfi.isValid(); ++mfi) {
          GridPtrDesc grid(leveldata, mfi);
          const Array4<CCTK_REAL> &vars = mfab->array(mfi);
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            CCTK_REAL *restrict const ptr = grid.ptr(vars, vi);
            grid.loop_all(groupdata.indextype, [&](const Loop::PointDesc &p) {
              ptr[p.idx] = 0.0 / 0.0;
            });
          }
        }
      }
      if (tl < prolongate_tl) {
        // Only interpolate if coarse grid data are valid
        bool all_invalid = true;
        for (int vi = 0; vi < groupdata.numvars; ++vi)
          all_invalid &= !coarsegroupdata.valid.at(tl).at(vi).valid_int &&
                         !coarsegroupdata.valid.at(tl).at(vi).valid_bnd &&
                         !groupdata.valid.at(tl).at(vi).valid_int &&
                         !groupdata.valid.at(tl).at(vi).valid_bnd;

        if (all_invalid) {
          valid.at(tl).valid_int = false;
          valid.at(tl).valid_bnd = false;
        } else {
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            const bool cond = coarsegroupdata.valid.at(tl).at(vi).valid_int &&
                              coarsegroupdata.valid.at(tl).at(vi).valid_bnd &&
                              groupdata.valid.at(tl).at(vi).valid_int &&
                              groupdata.valid.at(tl).at(vi).valid_bnd;
            if (!cond)
              CCTK_VERROR("Found invalid input data: RemakeLevel level %d, "
                          "variable %s%s: need everything defined, have coarse "
                          "%s, have current %s",
                          leveldata.level,
                          CCTK_FullVarName(groupdata.firstvarindex + vi),
                          string("_p", tl).c_str(),
                          string(coarsegroupdata.valid.at(tl).at(vi)).c_str(),
                          string(groupdata.valid.at(tl).at(vi)).c_str());
          }

          FillPatchTwoLevels(*mfab, 0.0, {&*coarsegroupdata.mfab.at(tl)}, {0.0},
                             {&*groupdata.mfab.at(tl)}, {0.0}, 0, 0,
                             groupdata.numvars, ghext->amrcore->Geom(level - 1),
                             ghext->amrcore->Geom(level), cphysbc, 0, fphysbc,
                             0, reffact, interpolator, bcs, 0);

          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            valid.at(tl).valid_int =
                coarsegroupdata.valid.at(tl).at(vi).valid_int &&
                coarsegroupdata.valid.at(tl).at(vi).valid_bnd &&
                groupdata.valid.at(tl).at(vi).valid_int &&
                groupdata.valid.at(tl).at(vi).valid_bnd;
            valid.at(tl).valid_bnd = false;
          }
        }
      }
      groupdata.mfab.at(tl) = move(mfab);
      groupdata.valid.at(tl) = move(valid);
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        poison_invalid(leveldata, groupdata, vi, tl);
    }
  }

  if (saved_cctkGH) {
    assert(current_level == -1);
    current_level = level;
    CCTK_Traverse(saved_cctkGH, "CCTK_BASEGRID");
    CCTK_Traverse(saved_cctkGH, "CCTK_POSTREGRID");
    current_level = -1;
  }
}

void CactusAmrCore::ClearLevel(int level) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("ClearLevel level %d", level);

  // assert(level == int(ghext->leveldata.size()) - 1);
  ghext->leveldata.resize(level);
}

////////////////////////////////////////////////////////////////////////////////

// Start driver
extern "C" int AMReX_Startup() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Startup");

  // Output a startup message
  string banner = "AMR driver provided by AMReX " + amrex::Version();
  int ierr = CCTK_RegisterBanner(banner.c_str());
  assert(!ierr);

  // Register a GH extension
  ghext_handle = CCTK_RegisterGHExtension("AMReX");
  assert(ghext_handle >= 0);
  int iret = CCTK_RegisterGHExtensionSetupGH(ghext_handle, SetupGH);
  assert(iret);
  iret = CCTK_RegisterGHExtensionInitGH(ghext_handle, InitGH);
  assert(iret);
  iret = CCTK_RegisterGHExtensionScheduleTraverseGH(ghext_handle,
                                                    ScheduleTraverseGH);
  assert(iret);

  CCTK_OverloadInitialise(Initialise);
  CCTK_OverloadEvolve(Evolve);
  CCTK_OverloadShutdown(Shutdown);
  CCTK_OverloadOutputGH(OutputGH);

  CCTK_OverloadMyProc(MyProc);
  CCTK_OverloadnProcs(nProcs);
  CCTK_OverloadExit(Exit);
  CCTK_OverloadAbort(Abort);
  CCTK_OverloadBarrier(Barrier);

  CCTK_OverloadSyncGroupsByDirI(SyncGroupsByDirI);

  return 0;
}

// Set up GH extension
void *SetupGH(tFleshConfig *fc, int convLevel, cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("SetupGH");

  assert(fc);
  assert(convLevel == 0);
  assert(cctkGH);

  // Initialize AMReX
  ParmParse pp;
  // Don't catch Unix signals. If signals are caught, we don't get
  // core files.
  pp.add("amrex.signal_handling", 0);
  // Throw exceptions for failing AMReX assertions. With exceptions,
  // we get core files.
  pp.add("amrex.throw_exception", 1);
  pamrex = amrex::Initialize(MPI_COMM_WORLD);

  // Create grid structure
  ghext = make_unique<GHExt>();

  return ghext.get();
}

// Initialize GH extension
int InitGH(cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("InitGH");

  assert(cctkGH);

  // Domain
  const RealBox domain({xmin, ymin, zmin}, {xmax, ymax, zmax});

  // Number of coarse grid cells
  const Vector<int> ncells{ncells_x, ncells_y, ncells_z};

  const int coord = -1; // undefined?

  // Refinement ratios
  const Vector<IntVect> reffacts; // empty

  // Periodicity
  const Array<int, dim> is_periodic{periodic_x, periodic_y, periodic_z};

  // Set blocking factors via parameter table since AmrMesh needs to
  // know them when its constructor is running, but there are no
  // constructor arguments for them
  ParmParse pp;
  pp.add("amr.blocking_factor_x", blocking_factor_x);
  pp.add("amr.blocking_factor_y", blocking_factor_y);
  pp.add("amr.blocking_factor_z", blocking_factor_z);
  pp.add("amr.max_grid_size_x", max_grid_size_x);
  pp.add("amr.max_grid_size_y", max_grid_size_y);
  pp.add("amr.max_grid_size_z", max_grid_size_z);
  pp.add("amr.grid_eff", grid_efficiency);

  ghext->amrcore = make_unique<CactusAmrCore>(
      domain, max_num_levels - 1, ncells, coord, reffacts, is_periodic);

  if (verbose) {
    int maxnumlevels = ghext->amrcore->maxLevel() + 1;
    for (int level = 0; level < maxnumlevels; ++level) {
      CCTK_VINFO("Geometry level %d:", level);
      cout << ghext->amrcore->Geom(level) << "\n";
    }
  }

  // CCTK_VINFO("BoxArray level %d:", level);
  // cout << ghext->amrcore->boxArray(level) << "\n";
  // CCTK_VINFO("DistributionMap level %d:", level);
  // cout << ghext->amrcore->DistributionMap(level) << "\n";

  return 0; // unused
}

// Traverse schedule
int ScheduleTraverseGH(cGH *restrict cctkGH, const char *where) {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("ScheduleTraverseGH [%d] %s", cctkGH->cctk_iteration, where);

  int ierr = CCTK_ScheduleTraverse(where, cctkGH, CallFunction);
  assert(!ierr);

  return 0; // unused
}

// Shut down driver
extern "C" int AMReX_Shutdown() {
  DECLARE_CCTK_PARAMETERS;
  if (verbose)
    CCTK_VINFO("Shutdown");

  // Should we really do this? Cactus's extension handling mechanism
  // becomes inconsistent once extensions have been unregistered.
  int iret = CCTK_UnregisterGHExtension("AMReX");
  assert(iret == 0);

  // Deallocate grid hierarchy
  ghext = nullptr;

  // Finalize AMReX
  amrex::Finalize(pamrex);
  pamrex = nullptr;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int MyProc(const cGH *restrict cctkGH) { return ParallelDescriptor::MyProc(); }

int nProcs(const cGH *restrict cctkGH) { return ParallelDescriptor::NProcs(); }

int Exit(cGH *cctkGH, int retval) {
  ParallelDescriptor::Abort(retval);
  return 0; // unreachable
}

int Abort(cGH *cctkGH, int retval) {
  ParallelDescriptor::Abort(retval);
  return 0; // unreachable
}

int Barrier(const cGH *restrict cctkGH) {
  ParallelDescriptor::Barrier();
  return 0;
}

} // namespace AMReX
