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

namespace {
// Convert a (direction, face) pair to an AMReX Orientation
Orientation orient(int d, int f) {
  return Orientation(d, Orientation::Side(f));
}
} // namespace

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
  // Don't regrid before Cactus is ready to
  if (level >= int(ghext->leveldata.size()))
    return;

  CCTK_VINFO("ErrorEst level %d", level);

  // // refine everywhere
  // tags.setVal(boxArray(level), TagBox::SET);

  // // refine centre
  // const Box &dom = Geom(level).Domain();
  // Box nbx;
  // for (int d = 0; d < dim; ++d) {
  //   int md = (dom.bigEnd(d) + dom.smallEnd(d) + 1) / 2;
  //   int rd = (dom.bigEnd(d) - dom.smallEnd(d) + 1) / 2;
  //   // mark one fewer cells; AMReX seems to add one cell
  //   nbx.setSmall(d, md - rd / (1 << (level + 1)) + 1);
  //   nbx.setBig(d, md + rd / (1 << (level + 1)) - 2);
  // }
  // cout << "EE nbx: " << nbx << "\n";
  // const BoxArray &ba = boxArray(level);
  // tags.setVal(intersect(ba, nbx), TagBox::SET);

  const int gi = CCTK_GroupIndex("AMReX::regrid_tag");
  assert(gi >= 0);
  const int vi = 0;
  const int tl = 0;

  auto &restrict leveldata = ghext->leveldata.at(level);
  auto &restrict groupdata = leveldata.groupdata.at(gi);
  MultiFab &mfab = *groupdata.mfab.at(tl);
  auto mfitinfo = MFItInfo().SetDynamic(true).EnableTiling({1024000, 16, 32});
#pragma omp parallel
  for (MFIter mfi(mfab, mfitinfo); mfi.isValid(); ++mfi) {
    const Box &fbx = mfi.fabbox();
    const Box &bx = mfi.tilebox();

    const Dim3 imin = lbound(bx);
    int ash[3], lsh[3];
    for (int d = 0; d < dim; ++d)
      ash[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)] + 1;
    for (int d = 0; d < dim; ++d)
      lsh[d] = bx[orient(d, 1)] - bx[orient(d, 0)] + 1;

    const Array4<CCTK_REAL> &ctagarr = groupdata.mfab.at(tl)->array(mfi);
    const CCTK_REAL *restrict ctags = ctagarr.ptr(imin.x, imin.y, imin.z, vi);
    const Array4<char> &atagarr = tags.array(mfi);
    char *restrict atags = atagarr.ptr(imin.x, imin.y, imin.z, vi);

    constexpr int di = 1;
    const int dj = di * ash[0];
    const int dk = dj * ash[1];
    for (int k = 0; k < lsh[2]; ++k) {
      for (int j = 0; j < lsh[1]; ++j) {
#pragma omp simd
        for (int i = 0; i < lsh[0]; ++i) {
          int idx = di * i + dj * j + dk * k;

          atags[idx] = ctags[idx] == 0.0 ? TagBox::CLEAR : TagBox::SET;
        }
      }
    }
  }
}

void CactusAmrCore::MakeNewLevelFromScratch(int lev, Real time,
                                            const BoxArray &ba,
                                            const DistributionMapping &dm) {
  CCTK_VINFO("MakeNewLevelFromScratch level %d", lev);
}

void CactusAmrCore::MakeNewLevelFromCoarse(int lev, Real time,
                                           const BoxArray &ba,
                                           const DistributionMapping &dm) {
  CCTK_VINFO("MakeNewLevelFromCoarse level %d", lev);
}

void CactusAmrCore::RemakeLevel(int lev, Real time, const BoxArray &ba,
                                const DistributionMapping &dm) {
  CCTK_VINFO("RemakeLevel level %d", lev);
}

void CactusAmrCore::ClearLevel(int lev) {
  CCTK_VINFO("ClearLevel level %d", lev);
}

////////////////////////////////////////////////////////////////////////////////

void SetupLevel(int level);

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

  // Number of coarse grid cells
  const Vector<int> ncells{ncells_x, ncells_y, ncells_z};

  const int coord = -1; // undefined?

  // Refinement ratios
  const Vector<IntVect> reffacts; // empty

  // Periodic in all directions
  const Array<int, dim> periodic{1, 1, 1};

  ghext->amrcore = make_unique<CactusAmrCore>(
      domain, max_num_levels - 1, ncells, coord, reffacts, periodic);
#warning "TODO: increase blocking factor"
  const int blocking_factor = 8;
  // const int blocking_factor = 1;
  ghext->amrcore->SetBlockingFactor(blocking_factor);
  const int max_grid_size = 32;
  ghext->amrcore->SetMaxGridSize(max_grid_size);

  int maxnumlevels = ghext->amrcore->maxLevel() + 1;
  for (int level = 0; level < maxnumlevels; ++level) {
    CCTK_VINFO("Geometry level %d:", level);
    cout << ghext->amrcore->Geom(level) << "\n";
  }

  // Create coarse grid
  const int level = 0;
  CCTK_REAL time = 0.0; // dummy time
  ghext->amrcore->MakeNewGrids(time);
  SetupLevel(level);

  // CCTK_VINFO("BoxArray level %d:", level);
  // cout << ghext->amrcore->boxArray(level) << "\n";
  // CCTK_VINFO("DistributionMap level %d:", level);
  // cout << ghext->amrcore->DistributionMap(level) << "\n";

  return 0; // unused
} // namespace AMReX

void CreateRefinedGrid(int level) {
  CCTK_VINFO("CreateRefinedGrid level %d", level);

  if (level > ghext->amrcore->maxLevel())
    return;

  // Create refined grid
  CCTK_REAL time = 0.0; // dummy time
  ghext->amrcore->regrid(0, time);
  int numlevels = ghext->amrcore->finestLevel() + 1;
  int maxnumlevels = ghext->amrcore->maxLevel() + 1;
  cout << "CRG numlevels=" << numlevels << "\n";
  cout << "CRG maxnumlevels=" << maxnumlevels << "\n";
  assert(numlevels >= 0 && numlevels <= maxnumlevels);
  assert(numlevels <= level + 1);

  if (numlevels == level + 1)
    SetupLevel(level);
}

void SetupLevel(int level) {
  DECLARE_CCTK_PARAMETERS;
  CCTK_VINFO("SetupLevel level %d", level);

  assert(level == int(ghext->leveldata.size()));
  ghext->leveldata.resize(level + 1);
  GHExt::LevelData &leveldata = ghext->leveldata.at(level);
  leveldata.level = level;

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
          ghext->amrcore->boxArray(leveldata.level),
          ghext->amrcore->DistributionMap(leveldata.level), groupdata.numvars,
          ghost_size);
    }
  }
}

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
