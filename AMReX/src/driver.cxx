#include <AMReX.hxx>
#include <schedule.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>

#include <AMReX.H>

#include <omp.h>
#include <mpi.h>

#include <string>
#include <type_traits>
#include <utility>

namespace AMReX {
using namespace amrex;
using namespace std;

// Global variables

int ghext_handle = -1;

amrex::AMReX *restrict pamrex = nullptr;
unique_ptr<GHExt> ghext;

vector<MFIter *> mfis;

// Registered functions

void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cctkGH);
int InitGH(cGH *cctkGH);
int ScheduleTraverseGH(cGH *cctkGH, const char *where);

int CallFunction(void *function, cFunctionData *attribute, void *data);

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
  CCTK_VINFO("InitGH");

  assert(cctkGH);

  // Define box array
  IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
  IntVect dom_hi(
      AMREX_D_DECL(ghext->ncells - 1, ghext->ncells - 1, ghext->ncells - 1));
  Box domain(dom_lo, dom_hi);
  ghext->ba.define(domain);

  // Break up box array into chunks no larger than max_grid_size along
  // each direction
  const int max_grid_size = 32;
  ghext->ba.maxSize(max_grid_size);

  // Define physical box
  RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                   {AMREX_D_DECL(1.0, 1.0, 1.0)});

  // Define geometry
  Vector<int> is_periodic(AMREX_SPACEDIM, 1); // periodic in all directions
  ghext->geom.define(domain, &real_box, CoordSys::cartesian,
                     is_periodic.data());

  CCTK_VINFO("Groups: %d", CCTK_NumGroups());
  CCTK_VINFO("Variables: %d", CCTK_NumVars());

  const int nvars = 4; // number of grid functions

  // Distributed boxes
  DistributionMapping dm(ghext->ba);

  // Allocate grid hierarchy
  ghext->mfab = MultiFab(ghext->ba, dm, nvars, ghext->nghostzones);

  setup_cctkGH(cctkGH);
  enter_global_mode(cctkGH);
  enter_level_mode(cctkGH);

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

} // namespace AMReX
