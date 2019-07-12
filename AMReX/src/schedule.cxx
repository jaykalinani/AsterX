#include <AMReX.hxx>
#include <schedule.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>

#include <AMReX.H>
#include <AMReX_Orientation.H>

#include <omp.h>
#include <mpi.h>

#include <string>
#include <type_traits>
#include <utility>

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

// Initialize cctkGH entries
void setup_cctkGH(cGH *restrict cctkGH) {
  // Grid function alignment
  // TODO: Check whether AMReX guarantees a particular alignment
  cctkGH->cctk_alignment = 1;
  cctkGH->cctk_alignment_offset = 0;

  // The refinement factor in time over the top level (coarsest) grid
  cctkGH->cctk_timefac = 1; // no subcycling

  // The convergence level (numbered from zero upwards)
  cctkGH->cctk_convlevel = 0; // no convergence tests

  // Initialize grid spacing
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_origin_space[d] = 0.0;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_delta_space[d] = 1.0;
  cctkGH->cctk_time = 0.0;
  cctkGH->cctk_delta_time = 0.0;
}

// Set cctkGH entries for global mode
void enter_global_mode(cGH *restrict cctkGH) {
  // The number of ghostzones in each direction
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_nghostzones[d] = ghext->nghostzones;
}
void leave_global_mode(cGH *restrict cctkGH) {
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_nghostzones[d] = undefined;
}

// Set cctkGH entries for local mode
void enter_level_mode(cGH *restrict cctkGH) {
  // Global shape
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_gsh[d] = ghext->ncells;

  // The refinement factor over the top level (coarsest) grid
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levfac[d] = 1; // TODO

  // Offset between this level's and the coarsest level's origin
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levoff[d] = 0; // TODO
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levoffdenom[d] = 0; // TODO
}
void leave_level_mode(cGH *restrict cctkGH) {
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_gsh[d] = undefined;
}

// Set cctkGH entries for local mode
void enter_local_mode(cGH *restrict cctkGH, const MFIter &mfi) {
  const Box &fbx = mfi.fabbox();   // allocated array
  const Box &vbx = mfi.validbox(); // interior region (without ghosts)
  const Box &bx = mfi.tilebox();   // current region (without ghosts)

  // Local shape
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lsh[d] = bx[orient(d, 1)] - bx[orient(d, 0)];

  // Allocated shape
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ash[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)];

  // Local extent
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lbnd[d] = bx[orient(d, 0)];
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ubnd[d] = bx[orient(d, 1)];

  // Boundaries
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      cctkGH->cctk_bbox[2 * d + f] = bx[orient(d, f)] == vbx[orient(d, f)];
}
void leave_local_mode(cGH *restrict cctkGH) {
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lsh[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ash[d] = undefined;
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
    mfis.resize(omp_get_max_threads(), nullptr);
#pragma omp parallel
    for (MFIter mfi(ghext->mfab, MFItInfo().SetDynamic(true).EnableTiling(
                                     {1024000, 16, 32}));
         mfi.isValid(); ++mfi) {
      mfis.at(omp_get_thread_num()) = &mfi;
      enter_local_mode(cctkGH, mfi);
      CCTK_CallFunction(function, attribute, data);
    }
    mfis.clear();
    leave_local_mode(cctkGH);
    break;
  }
  case mode_t::level:
  case mode_t::global:
  case mode_t::meta: {
    // Call function once
    // Note: meta mode scheduling must continue to work even after we
    // shut down ourselves!
    CCTK_CallFunction(function, attribute, data);
    break;
  }
  default:
    assert(0);
  }

  int didsync = 1;
  return didsync;
}

} // namespace AMReX
