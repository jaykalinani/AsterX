#include <AMReX.hxx>
#include <schedule.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>
#include <cctki_GHExtensions.h>
#include <cctki_ScheduleBindings.h>
#include <cctki_WarnLevel.h>

#include <AMReX.H>
#include <AMReX_Orientation.H>

#include <omp.h>
#include <mpi.h>

#include <sys/time.h>

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace AMReX {
using namespace amrex;
using namespace std;

// Value for undefined cctkGH entries
// Note: Don't use a negative value, which tends to leave bugs undetected. Large
// positive values often lead to segfault, exposing bugs.
constexpr int undefined = 666;

////////////////////////////////////////////////////////////////////////////////

// Convert a (direction, face) pair to an AMReX Orientation
Orientation orient(int d, int f) {
  return Orientation(d, Orientation::Side(f));
}

vector<cGH> thread_local_cctkGH;

// Create a new cGH, copying those data that are set by the flesh, and
// allocating space for these data that are set per thread by the driver
void clone_cctkGH(cGH *restrict cctkGH, const cGH *restrict sourceGH) {
  // Copy all fields by default
  *cctkGH = *sourceGH;
  // Allocate most fields anew
  cctkGH->cctk_gsh = new int[dim];
  cctkGH->cctk_lsh = new int[dim];
  cctkGH->cctk_lbnd = new int[dim];
  cctkGH->cctk_ubnd = new int[dim];
  cctkGH->cctk_ash = new int[dim];
  cctkGH->cctk_to = new int[dim];
  cctkGH->cctk_from = new int[dim];
  cctkGH->cctk_delta_space = new CCTK_REAL[dim];
  cctkGH->cctk_origin_space = new CCTK_REAL[dim];
  cctkGH->cctk_bbox = new int[2 * dim];
  cctkGH->cctk_levfac = new int[dim];
  cctkGH->cctk_levoff = new int[dim];
  cctkGH->cctk_levoffdenom = new int[dim];
  cctkGH->cctk_nghostzones = new int[dim];
  const int numvars = CCTK_NumVars();
  cctkGH->data = new void **[numvars];
  for (int vi = 0; vi < numvars; ++vi)
    cctkGH->data[vi] = new void *[CCTK_DeclaredTimeLevelsVI(vi)];
}

// Initialize cctkGH entries
void setup_cctkGH(cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // Grid function alignment
  // TODO: Check whether AMReX guarantees a particular alignment
  cctkGH->cctk_alignment = 1;
  cctkGH->cctk_alignment_offset = 0;

  // The refinement factor in time over the top level (coarsest) grid
  cctkGH->cctk_timefac = 1; // no subcycling

  // The convergence level (numbered from zero upwards)
  cctkGH->cctk_convlevel = 0; // no convergence tests

  // Initialize grid spacing
  int ncells[dim] = {ncells_x, ncells_y, ncells_z};
  CCTK_REAL x0[dim] = {xmin, ymin, zmin};
  CCTK_REAL x1[dim] = {xmax, ymax, zmax};
  CCTK_REAL dx[dim];
  for (int d = 0; d < dim; ++d)
    dx[d] = (x1[d] - x0[d]) / ncells[d];
  CCTK_REAL mindx = 1.0 / 0.0;
  for (int d = 0; d < dim; ++d)
    mindx = fmin(mindx, dx[d]);

  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_origin_space[d] = x0[d] + 0.5 * dx[d];
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_delta_space[d] = dx[d];

  // Initialize time stepping
  cctkGH->cctk_time = 0.0;
  cctkGH->cctk_delta_time = dtfac * mindx;
}

// Update fields that carry state and change over time
void update_cctkGH(cGH *restrict cctkGH, const cGH *restrict sourceGH) {
  cctkGH->cctk_iteration = sourceGH->cctk_iteration;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_origin_space[d] = sourceGH->cctk_origin_space[d];
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_delta_space[d] = sourceGH->cctk_delta_space[d];
  cctkGH->cctk_time = sourceGH->cctk_time;
  cctkGH->cctk_delta_time = sourceGH->cctk_delta_time;
}

// Set cctkGH entries for global mode
void enter_global_mode(cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  // The number of ghostzones in each direction
  // TODO: Get this from mfab
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_nghostzones[d] = ghost_size;
}
void leave_global_mode(cGH *restrict cctkGH) {
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_nghostzones[d] = undefined;
}

// Set cctkGH entries for local mode
void enter_level_mode(cGH *restrict cctkGH,
                      const GHExt::Level &restrict level) {
  DECLARE_CCTK_PARAMETERS;

  // Global shape
  // TODO: Get this from mfab
#warning "TODO"
  int ncells[dim] = {ncells_x, ncells_y, ncells_z};
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_gsh[d] = ncells[d];

    // The refinement factor over the top level (coarsest) grid
#warning "TODO"
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levfac[d] = 1; // TODO

    // Offset between this level's and the coarsest level's origin
#warning "TODO"
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levoff[d] = 0; // TODO
#warning "TODO"
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levoffdenom[d] = 0; // TODO
}
void leave_level_mode(cGH *restrict cctkGH,
                      const GHExt::Level &restrict level) {
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_gsh[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levfac[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levoff[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levoffdenom[d] = 0;
}

// Set cctkGH entries for local mode
void enter_local_mode(cGH *restrict cctkGH, GHExt::Level &restrict level,
                      const MFIter &mfi) {
  const Box &fbx = mfi.fabbox(); // allocated array
  // const Box &vbx = mfi.validbox();     // interior region (without ghosts)
  const Box &bx = mfi.tilebox();       // current region (without ghosts)
  const Box &gbx = mfi.growntilebox(); // current region (with ghosts)

  // Local shape
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lsh[d] = bx[orient(d, 1)] - bx[orient(d, 0)] + 1;

  // Allocated shape
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ash[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)] + 1;

  // Local extent
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lbnd[d] = bx[orient(d, 0)];
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ubnd[d] = bx[orient(d, 1)];

  // Boundaries
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      cctkGH->cctk_bbox[2 * d + f] = bx[orient(d, f)] != gbx[orient(d, f)];

  // Grid function pointers
  const Dim3 imin = lbound(bx);
  for (auto &restrict groupdata : level.groupdata) {
    const Array4<CCTK_REAL> &vars = groupdata.mfab.array(mfi);
    for (int tl = 0; tl < groupdata.numtimelevels; ++tl)
      for (int n = 0; n < groupdata.numvars; ++n)
        cctkGH->data[groupdata.firstvarindex + n][tl] =
            vars.ptr(imin.x, imin.y, imin.z, tl * groupdata.numvars + n);
  }
}
void leave_local_mode(cGH *restrict cctkGH, const GHExt::Level &restrict level,
                      const MFIter &mfi) {
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lsh[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ash[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lbnd[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ubnd[d] = undefined;
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      cctkGH->cctk_bbox[2 * d + f] = undefined;
  for (auto &restrict groupdata : level.groupdata) {
    for (int tl = 0; tl < groupdata.numtimelevels; ++tl)
      for (int n = 0; n < groupdata.numvars; ++n)
        cctkGH->data[groupdata.firstvarindex + n][tl] = nullptr;
  }
}

enum class mode_t { unknown, local, level, global, meta };

mode_t decode_mode(const cFunctionData *restrict attribute) {
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

// Schedule initialisation
int Initialise(tFleshConfig *config) {
  cGH *restrict const cctkGH = CCTK_SetupGH(config, 0);
  CCTKi_AddGH(config, 0, cctkGH);

  // Initialise iteration and time
  cctkGH->cctk_iteration = 0;
  cctkGH->cctk_time = *static_cast<const CCTK_REAL *>(
      CCTK_ParameterGet("cctk_initial_time", "Cactus", nullptr));

  // Initialise schedule
  CCTKi_ScheduleGHInit(cctkGH);

  // Initialise all grid extensions
  CCTKi_InitGHExtensions(cctkGH);

  // Set up cctkGH
  setup_cctkGH(cctkGH);
  enter_global_mode(cctkGH);
  // enter_level_mode(cctkGH);
  int max_threads = omp_get_max_threads();
  thread_local_cctkGH.resize(max_threads);
  for (int n = 0; n < max_threads; ++n) {
    cGH *restrict threadGH = &thread_local_cctkGH.at(n);
    clone_cctkGH(threadGH, cctkGH);
    setup_cctkGH(threadGH);
    enter_global_mode(threadGH);
    // enter_level_mode(threadGH);
  }

  // Output domain information
  if (CCTK_MyProc(nullptr) == 0) {
    int ncells[dim];
    for (int d = 0; d < dim; ++d)
      ncells[d] = cctkGH->cctk_gsh[d];
    CCTK_REAL x0[dim], x1[dim], dx[dim];
    for (int d = 0; d < dim; ++d)
      dx[d] = cctkGH->cctk_delta_space[d];
    for (int d = 0; d < dim; ++d)
      x0[d] = cctkGH->cctk_origin_space[d] - 0.5 * dx[d];
    for (int d = 0; d < dim; ++d)
      x1[d] = x0[d] + ncells[d] * dx[d];
    CCTK_VINFO("Grid extent:");
    CCTK_VINFO("  gsh=[%d,%d,%d]", ncells[0], ncells[1], ncells[2]);
    CCTK_VINFO("Domain extent:");
    CCTK_VINFO("  xmin=[%g,%g,%g]", x0[0], x0[1], x0[2]);
    CCTK_VINFO("  xmax=[%g,%g,%g]", x1[0], x1[1], x1[2]);
    CCTK_VINFO("  dx=[%g,%g,%g]", dx[0], dx[1], dx[2]);
    CCTK_VINFO("Time stepping:");
    CCTK_VINFO("  t0=%g", cctkGH->cctk_time);
    CCTK_VINFO("  dt=%g", cctkGH->cctk_delta_time);
  }

  CCTK_Traverse(cctkGH, "CCTK_WRAGH");
  CCTK_Traverse(cctkGH, "CCTK_PARAMCHECK");
  CCTKi_FinaliseParamWarn();

  CCTK_Traverse(cctkGH, "CCTK_BASEGRID");

  const char *recovery_mode = *static_cast<const char *const *>(
      CCTK_ParameterGet("recovery_mode", "Cactus", nullptr));
  if (!config->recovered || !CCTK_Equals(recovery_mode, "strict")) {
    // Set up initial conditions
    CCTK_Traverse(cctkGH, "CCTK_INITIAL");
    CCTK_Traverse(cctkGH, "CCTK_POSTINITIAL");
    CCTK_Traverse(cctkGH, "CCTK_POSTPOSTINITIAL");
    CCTK_Traverse(cctkGH, "CCTK_POSTSTEP");
  }

  if (config->recovered) {
    // Recover
    CCTK_Traverse(cctkGH, "CCTK_RECOVER_VARIABLES");
    CCTK_Traverse(cctkGH, "CCTK_POST_RECOVER_VARIABLES");
  }

  // Checkpoint, analysis, output
  CCTK_Traverse(cctkGH, "CCTK_CPINITIAL");
  CCTK_Traverse(cctkGH, "CCTK_ANALYSIS");
  CCTK_OutputGH(cctkGH);

  return 0;
}

bool EvolutionIsDone(cGH *restrict const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  static timeval starttime = {0, 0};
  // On the first time through, get the start time
  if (starttime.tv_sec == 0 && starttime.tv_usec == 0)
    gettimeofday(&starttime, nullptr);

  if (terminate_next || CCTK_TerminationReached(cctkGH))
    return true;

  if (CCTK_Equals(terminate, "never"))
    return false;

  bool max_iteration_reached = cctkGH->cctk_iteration >= cctk_itlast;

  bool max_simulation_time_reached = cctk_initial_time < cctk_final_time
                                         ? cctkGH->cctk_time >= cctk_final_time
                                         : cctkGH->cctk_time <= cctk_final_time;

  // Get the elapsed runtime in minutes and compare with max_runtime
  timeval currenttime;
  gettimeofday(&currenttime, NULL);
  bool max_runtime_reached =
      CCTK_REAL(currenttime.tv_sec - starttime.tv_sec) / 60 >= max_runtime;

  if (CCTK_Equals(terminate, "iteration"))
    return max_iteration_reached;
  if (CCTK_Equals(terminate, "time"))
    return max_simulation_time_reached;
  if (CCTK_Equals(terminate, "runtime"))
    return max_runtime_reached;
  if (CCTK_Equals(terminate, "any"))
    return max_iteration_reached || max_simulation_time_reached ||
           max_runtime_reached;
  if (CCTK_Equals(terminate, "all"))
    return max_iteration_reached && max_simulation_time_reached &&
           max_runtime_reached;
  if (CCTK_Equals(terminate, "either"))
    return max_iteration_reached || max_simulation_time_reached;
  if (CCTK_Equals(terminate, "both"))
    return max_iteration_reached && max_simulation_time_reached;

  assert(0);
}

void CycleTimelevels(cGH *restrict const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  cctkGH->cctk_iteration += 1;
  cctkGH->cctk_time += cctkGH->cctk_delta_time;

  // TODO: Get ghost_size from mfab
  for (auto &restrict level : ghext->levels)
    for (auto &restrict groupdata : level.groupdata) {
      for (int tl = groupdata.numtimelevels - 1; tl > 0; --tl)
        MultiFab::Copy(groupdata.mfab, groupdata.mfab,
                       (tl - 1) * groupdata.numvars, tl * groupdata.numvars,
                       groupdata.numvars, ghost_size);
    }
}

// Schedule evolution
int Evolve(tFleshConfig *config) {
  assert(config);
  cGH *restrict const cctkGH = config->GH[0];
  assert(cctkGH);

  while (!EvolutionIsDone(cctkGH)) {
    CycleTimelevels(cctkGH);

    CCTK_Traverse(cctkGH, "CCTK_PRESTEP");
    CCTK_Traverse(cctkGH, "CCTK_EVOL");
    CCTK_Traverse(cctkGH, "CCTK_POSTSTEP");

    CCTK_Traverse(cctkGH, "CCTK_CHECKPOINT");
    CCTK_Traverse(cctkGH, "CCTK_ANALYSIS");
    CCTK_OutputGH(cctkGH);
  }

  return 0;
}

// Schedule shutdown
int Shutdown(tFleshConfig *config) {
  assert(config);
  cGH *restrict const cctkGH = config->GH[0];
  assert(cctkGH);

  CCTK_Traverse(cctkGH, "CCTK_TERMINATE");
  CCTK_Traverse(cctkGH, "CCTK_SHUTDOWN");

  return 0;
}

// Call a scheduled function
int CallFunction(void *function, cFunctionData *restrict attribute,
                 void *data) {
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
#pragma omp parallel
    {
      cGH *restrict threadGH = &thread_local_cctkGH.at(omp_get_thread_num());
      update_cctkGH(threadGH, cctkGH);
      for (auto &restrict level : ghext->levels) {
        MultiFab &mfab = level.groupdata.at(0).mfab;
        enter_level_mode(threadGH, level);
        auto mfitinfo =
            MFItInfo().SetDynamic(true).EnableTiling({1024000, 16, 32});
        for (MFIter mfi(mfab, mfitinfo); mfi.isValid(); ++mfi) {
          enter_local_mode(threadGH, level, mfi);
          CCTK_CallFunction(function, attribute, threadGH);
          leave_local_mode(threadGH, level, mfi);
        }
        leave_level_mode(threadGH, level);
      }
    }
    break;
  }
  case mode_t::meta:
  case mode_t::global:
  case mode_t::level: {
    // Call function just once
    // Note: meta mode scheduling must continue to work even after we
    // shut down ourselves!
    CCTK_CallFunction(function, attribute, cctkGH);
    break;
  }

  default:
    assert(0);
  }

  int didsync = 0;
  return didsync;
}

int SyncGroupsByDirI(const cGH *restrict cctkGH, int numgroups,
                     const int *groups, const int *directions) {
  assert(cctkGH);
  assert(numgroups >= 0);
  assert(groups);

  // TODO: Synchronize only current time level. This probably requires
  // setting up different mfabs for each time level.
  for (auto &restrict level : ghext->levels)
    for (int n = 0; n < numgroups; ++n) {
      int gi = groups[n];
      auto &restrict groupdata = level.groupdata.at(gi);
      // we always sync all directions
      groupdata.mfab.FillBoundary(level.geom.periodicity());
    }

  return numgroups; // number of groups synchronized
}

} // namespace AMReX
