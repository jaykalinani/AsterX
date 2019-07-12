#include <AMReX.hxx>

#include <AMReX.H>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>

#include <omp.h>
#include <mpi.h>

#include <string>
#include <type_traits>
#include <utility>

namespace AMReX {
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
void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cctkGH) {
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
int InitGH(cGH *cctkGH) {
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

  return 0; // unused
}

// Traverse schedule
int ScheduleTraverseGH(cGH *cctkGH, const char *where) {
  CCTK_VINFO("ScheduleTraverseGH [%d] %s", cctkGH->cctk_iteration, where);

  int ierr = CCTK_ScheduleTraverse(where, cctkGH, CallFunction);
  assert(!ierr);

  return 0; // unused
}

namespace {

enum class mode_t { unknown, local, level, global, meta };

mode_t decode_mode(const cFunctionData *attribute) {
  bool local_mode = attribute->local;
  bool level_mode = attribute->level;
  bool global_mode = attribute->global;
  bool meta_mode = attribute->meta;
  assert(int(local_mode) + int(level_mode) + int(global_mode) +
             int(meta_mode) <=
         1);
  if (attribute->local)
    return mode_t::local;
  if (attribute->level)
    return mode_t::level;
  if (attribute->global)
    return mode_t::global;
  if (attribute->meta)
    return mode_t::meta;
  return mode_t::local; // default
}
} // namespace

// Call a scheduled function
int CallFunction(void *function, cFunctionData *attribute, void *data) {
  assert(function);
  assert(attribute);
  assert(data);

  cGH *restrict const cctkGH = static_cast<cGH *>(data);

  CCTK_VINFO("CallFunction [%d] %s: %s::%s", cctkGH->cctk_iteration,
             attribute->where, attribute->thorn, attribute->routine);

  const mode_t mode = decode_mode(attribute);
  switch (mode) {
  case mode_t::local: {
    // Call function once per tile
    mfis.resize(omp_get_max_threads(), nullptr);
#pragma omp parallel
    for (MFIter mfi(ghext->mfab, MFItInfo().SetDynamic(true).EnableTiling(
                                     {1024000, 16, 32}));
         mfi.isValid(); ++mfi) {
      mfis.at(omp_get_thread_num()) = &mfi;
      CCTK_CallFunction(function, attribute, data);
    }
    mfis.clear();
    break;
  }
  case mode_t::level:
  case mode_t::global:
  case mode_t::meta: {
    // Call function once
    CCTK_CallFunction(function, attribute, data);
    break;
  }
  default:
    assert(0);
  }

  int didsync = 1;
  return didsync;
}

// Shut down driver
extern "C" int AMReX_Shutdown() {
  CCTK_VINFO("Shutdown");

  int iret = CCTK_UnregisterGHExtension("AMReX");
  CCTK_VINFO("iret=%d", iret);
  assert(iret == 0);

  // Deallocate grid hierarchy
  ghext = nullptr;

  // Finalize AMReX
  amrex::Finalize(pamrex);
  pamrex = nullptr;

  return 0;
}

} // namespace AMReX
