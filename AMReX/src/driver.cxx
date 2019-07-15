#include <AMReX.hxx>
#include <io.hxx>
#include <schedule.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>

#include <omp.h>
#include <mpi.h>

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

  // Define box array
  IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
  IntVect dom_hi(AMREX_D_DECL(ncells_x - 1, ncells_y - 1, ncells_z - 1));
  Box domain(dom_lo, dom_hi);
  ghext->ba.define(domain);

  // Break up box array into chunks no larger than max_grid_size along
  // each direction
  const int max_grid_size = 32;
  ghext->ba.maxSize(max_grid_size);

  // Define physical box
  RealBox real_box({AMREX_D_DECL(xmin, ymin, zmin)},
                   {AMREX_D_DECL(xmax, ymax, zmax)});

  // Define geometry
  Vector<int> is_periodic(AMREX_SPACEDIM, 1); // periodic in all directions
  ghext->geom.define(domain, &real_box, CoordSys::cartesian,
                     is_periodic.data());

  // Distributed boxes
  DistributionMapping dm(ghext->ba);

  const int numgroups = CCTK_NumGroups();
  ghext->groupdata.resize(numgroups);
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);
    assert(group.grouptype == CCTK_GF);
    assert(group.vartype == CCTK_VARIABLE_REAL);
    assert(group.disttype == CCTK_DISTRIB_DEFAULT);
    assert(group.dim == dim);

    GHExt::GroupData &groupdata = ghext->groupdata.at(gi);
    groupdata.firstvarindex = CCTK_FirstVarIndexI(gi);
    groupdata.numvars = group.numvars;
    groupdata.numtimelevels = group.numtimelevels;

    // Allocate grid hierarchy
    groupdata.mfab = MultiFab(
        ghext->ba, dm, groupdata.numvars * groupdata.numtimelevels, ghost_size);
  }

  return 0; // unused
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
