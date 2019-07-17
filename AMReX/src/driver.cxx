#include <driver.hxx>
#include <io.hxx>
#include <schedule.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>
#include <AMReX_BCRec.H>

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

CactusAmrMesh::CactusAmrMesh() {}
CactusAmrMesh::CactusAmrMesh(const RealBox *rb, int max_level_in,
                             const Vector<int> &n_cell_in, int coord,
                             Vector<IntVect> ref_ratios, const int *is_per)
    : AmrMesh(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per) {}
CactusAmrMesh::CactusAmrMesh(const RealBox &rb, int max_level_in,
                             const Vector<int> &n_cell_in, int coord,
                             Vector<IntVect> const &ref_ratios,
                             Array<int, AMREX_SPACEDIM> const &is_per)
    : AmrMesh(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per) {}

CactusAmrMesh::~CactusAmrMesh() {}

void CactusAmrMesh::ErrorEst(int lev, TagBoxArray &tags, Real time, int ngrow) {
  // // refine everywhere
  // tags.setVal(boxArray(lev), TagBox::SET);

  // refine centre
  const BoxArray &ba = boxArray(lev);
  const Box &bx = ba.minimalBox();
  Box nbx;
  for (int d = 0; d < dim; ++d) {
    int sz = bx.bigEnd(d) - bx.smallEnd(d) + 1;
    nbx.setSmall(d, bx.smallEnd(d) + sz / 4 + 1);
    nbx.setBig(d, bx.bigEnd(d) - sz / 4 - 1);
  }
  cout << "nbx: " << nbx << "\n";
  tags.setVal(intersect(ba, nbx), TagBox::SET);
}

////////////////////////////////////////////////////////////////////////////////

// Start driver
extern "C" int AMReX_Startup() {
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
  CCTK_VINFO("SetupGH");

  assert(fc);
  assert(convLevel == 0);
  assert(cctkGH);

  // Initialize AMReX
  pamrex = amrex::Initialize(MPI_COMM_WORLD);

  // Create grid structure
  ghext = make_unique<GHExt>();

  return ghext.get();
}

// Initialize GH extension
int InitGH(cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;
  CCTK_VINFO("InitGH");

  assert(cctkGH);

  // Domain
  const RealBox domain({xmin, ymin, zmin}, {xmax, ymax, zmax});

  // Maximum number of levels
  const int maxnumlevels = 2;

  // Number of coarse grid cells
  const Vector<int> ncells{ncells_x, ncells_y, ncells_z};

  const int coord = -1; // undefined?

  // Refinement ratios
  const Vector<IntVect> reffacts; // empty

  // Periodic in all directions
  const Array<int, dim> periodic{1, 1, 1};

  ghext->amrmesh = make_unique<CactusAmrMesh>(domain, maxnumlevels - 1, ncells,
                                              coord, reffacts, periodic);
#warning "TODO: increase blocking factor"
  // const int blocking_factor = 8;
  const int blocking_factor = 1;
  ghext->amrmesh->SetBlockingFactor(blocking_factor);
  const int max_grid_size = 32;
  ghext->amrmesh->SetMaxGridSize(max_grid_size);

  for (int level = 0; level < maxnumlevels; ++level) {
    CCTK_VINFO("Geometry level %d:", level);
    cout << ghext->amrmesh->Geom(level) << "\n";
  }

  const int numlevels = maxnumlevels;
  ghext->leveldata.resize(numlevels);
  for (int level = 0; level < numlevels; ++level) {
    GHExt::LevelData &leveldata = ghext->leveldata.at(level);
    leveldata.level = level;

    CCTK_REAL time = 0.0; // dummy time

    // Create grid structure
    if (level == 0) {
      // Create coarse grid
      ghext->amrmesh->MakeNewGrids(time);
    } else {
      // Create refined grid
      int new_finest = -999;
      Vector<BoxArray> new_grids;
      ghext->amrmesh->MakeNewGrids(0, time, new_finest, new_grids);
      assert(new_finest == level);
    }

    // CCTK_VINFO("BoxArray level %d:", level);
    // cout << ghext->amrmesh->boxArray(level) << "\n";
    // CCTK_VINFO("DistributionMap level %d:", level);
    // cout << ghext->amrmesh->DistributionMap(level) << "\n";

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
      groupdata.firstvarindex = CCTK_FirstVarIndexI(gi);
      groupdata.numvars = group.numvars;

      // Allocate grid hierarchies
      groupdata.mfab.resize(group.numtimelevels);
      for (int tl = 0; tl < int(groupdata.mfab.size()); ++tl) {
        groupdata.mfab.at(tl) = make_unique<MultiFab>(
            ghext->amrmesh->boxArray(leveldata.level),
            ghext->amrmesh->DistributionMap(leveldata.level), groupdata.numvars,
            ghost_size);
      }
    }

  } // for level

  return 0; // unused
} // namespace AMReX

// Traverse schedule
int ScheduleTraverseGH(cGH *restrict cctkGH, const char *where) {
  CCTK_VINFO("ScheduleTraverseGH [%d] %s", cctkGH->cctk_iteration, where);

  int ierr = CCTK_ScheduleTraverse(where, cctkGH, CallFunction);
  assert(!ierr);

  return 0; // unused
}

// Shut down driver
extern "C" int AMReX_Shutdown() {
  CCTK_VINFO("Shutdown");

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
