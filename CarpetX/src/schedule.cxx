#include "driver.hxx"
#include "io.hxx"
#include "loop.hxx"
#include "schedule.hxx"
#include "timer.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cctk_Schedule.h>
#include <cctki_GHExtensions.h>
#include <cctki_ScheduleBindings.h>
#include <cctki_WarnLevel.h>
#include <util_Table.h>

#include <AMReX.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Orientation.H>
#include <AMReX_PhysBCFunct.H>

#ifdef _OPENMP
#include <omp.h>
#else
extern "C" {
static inline int omp_get_max_threads(void) { return 1; }
static inline int omp_get_num_threads(void) { return 1; }
static inline int omp_get_thread_num(void) { return 0; }
static inline int omp_in_parallel(void) { return 0; }
}
#endif
#include <mpi.h>
#include <zlib.h>

#include <sys/time.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace CarpetX {
using namespace std;

// Value for undefined cctkGH entries
// Note: Don't use a negative value, which tends to leave bugs undetected. Large
// positive values often lead to segfault, exposing bugs.
constexpr int undefined = 666;

// Tile boxes, should be part of cGH
struct TileBox {
  array<int, dim> tile_min;
  array<int, dim> tile_max;
};

struct thread_local_info_t {
  // TODO: store only amrex::MFIter here; recalculate other things from it
  cGH cctkGH;
  TileBox tilebox;
  const MFPointer *restrict mfpointer;
  unsigned char padding[128]; // Prevent false sharing
};

vector<unique_ptr<thread_local_info_t> > thread_local_info;
vector<unique_ptr<thread_local_info_t> > saved_thread_local_info;

// Used to pass cctkGH into the callbacks for AMReX's regridding functions
cGH *saved_cctkGH = nullptr;

// Used to pass active levels from AMReX's regridding functions
optional<active_levels_t> active_levels;

void Reflux(int level);
void Restrict(int level, const vector<int> &groups);
void Restrict(int level);

namespace {
// Convert a (direction, face) pair to an AMReX Orientation
amrex::Orientation orient(int d, int f) {
  return amrex::Orientation(d, amrex::Orientation::Side(f));
}
int GroupStorageCrease(const cGH *cctkGH, int n_groups, const int *groups,
                       const int *requested_tls, int *status, const bool inc);
} // namespace

////////////////////////////////////////////////////////////////////////////////

// Set up a GridDesc
} // namespace CarpetX
namespace Loop {

GridDescBase::GridDescBase() {}

GridDescBase::GridDescBase(const cGH *restrict cctkGH) {
  for (int d = 0; d < dim; ++d)
    gsh[d] = cctkGH->cctk_gsh[d];
  for (int d = 0; d < dim; ++d) {
    lbnd[d] = cctkGH->cctk_lbnd[d];
    ubnd[d] = cctkGH->cctk_ubnd[d];
  }
  for (int d = 0; d < dim; ++d)
    lsh[d] = cctkGH->cctk_lsh[d];
  for (int d = 0; d < dim; ++d)
    ash[d] = cctkGH->cctk_ash[d];
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      bbox[2 * d + f] = cctkGH->cctk_bbox[2 * d + f];
  for (int d = 0; d < dim; ++d)
    nghostzones[d] = cctkGH->cctk_nghostzones[d];

  // Check whether we are in local mode
  assert(cctkGH->cctk_bbox[0] != CarpetX::undefined);
  int thread_num = omp_get_thread_num();
  if (omp_in_parallel()) {
    const cGH *restrict threadGH =
        &CarpetX::thread_local_info.at(thread_num)->cctkGH;
    // Check whether this is the correct cGH structure
    assert(cctkGH == threadGH);
  }

  const CarpetX::TileBox &restrict tilebox =
      CarpetX::thread_local_info.at(thread_num)->tilebox;
  for (int d = 0; d < dim; ++d) {
    tmin[d] = tilebox.tile_min[d];
    tmax[d] = tilebox.tile_max[d];
  }

  for (int d = 0; d < dim; ++d) {
    dx[d] = cctkGH->cctk_delta_space[d] / cctkGH->cctk_levfac[d];
    x0[d] = cctkGH->cctk_origin_space[d] +
            dx[d] * cctkGH->cctk_levoff[d] / cctkGH->cctk_levoffdenom[d];
  }
}

} // namespace Loop
namespace CarpetX {

GridDesc::GridDesc(const GHExt::LevelData &leveldata, const MFPointer &mfp) {
  DECLARE_CCTK_PARAMETERS;

  const amrex::Box &fbx = mfp.fabbox();   // allocated array
  const amrex::Box &vbx = mfp.validbox(); // interior region (without ghosts)
  const amrex::Box &gbx = mfp.growntilebox(); // current region (with ghosts)
  const amrex::Box &domain = ghext->amrcore->Geom(leveldata.level).Domain();

  // The number of ghostzones in each direction
  for (int d = 0; d < dim; ++d)
    nghostzones[d] = mfp.nGrowVect()[d];

  // Global shape
  for (int d = 0; d < dim; ++d)
    gsh[d] =
        domain[orient(d, 1)] + 1 - domain[orient(d, 0)] + 2 * nghostzones[d];

  // Local shape
  for (int d = 0; d < dim; ++d)
    lsh[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)] + 1;

  // Allocated shape
  for (int d = 0; d < dim; ++d)
    ash[d] = fbx[orient(d, 1)] - fbx[orient(d, 0)] + 1;

  // Local extent
  for (int d = 0; d < dim; ++d) {
    lbnd[d] = fbx[orient(d, 0)] + nghostzones[d];
    ubnd[d] = fbx[orient(d, 1)] + nghostzones[d];
  }

  // Boundaries
  const array<array<bool, 3>, 2> is_symmetry{{
      {{
          periodic || periodic_x || reflection_x,
          periodic || periodic_y || reflection_y,
          periodic || periodic_z || reflection_z,
      }},
      {{
          periodic || periodic_x || reflection_upper_x,
          periodic || periodic_y || reflection_upper_y,
          periodic || periodic_z || reflection_upper_z,
      }},
  }};
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      bbox[2 * d + f] =
          vbx[orient(d, f)] == domain[orient(d, f)] && !is_symmetry[f][d];

  // Thread tile box
  for (int d = 0; d < dim; ++d) {
    tmin[d] = gbx[orient(d, 0)] - fbx[orient(d, 0)];
    tmax[d] = gbx[orient(d, 1)] + 1 - fbx[orient(d, 0)];
  }

  const amrex::Geometry &geom = ghext->amrcore->Geom(0);
  const CCTK_REAL *restrict const global_x0 = geom.ProbLo();
  const CCTK_REAL *restrict const global_dx = geom.CellSize();
  for (int d = 0; d < dim; ++d) {
    const int levfac = 1 << leveldata.level;
    const int levoff = 1 - 2 * nghostzones[d];
    const int levoffdenom = 2;
    const CCTK_REAL origin_space =
        global_x0[d] + (1 - nghostzones[d] + CCTK_REAL(1) / 2) * global_dx[d];
    const CCTK_REAL delta_space = global_dx[d];
    dx[d] = delta_space / levfac;
    x0[d] = origin_space + dx[d] * levoff / levoffdenom;
  }

  // Check constraints
  for (int d = 0; d < dim; ++d) {
    // Domain size
    assert(gsh[d] >= 0);

    // Local size
    assert(lbnd[d] >= 0);
    assert(lsh[d] >= 0);
    assert(lbnd[d] + lsh[d] <= gsh[d]);
    assert(ubnd[d] == lbnd[d] + lsh[d] - 1);

    // Internal representation
    assert(ash[d] >= 0);
    assert(ash[d] >= lsh[d]);

    // Ghost zones
    assert(nghostzones[d] >= 0);
    assert(2 * nghostzones[d] <= lsh[d]);

    // Tiles
    assert(tmin[d] >= 0);
    assert(tmax[d] >= tmin[d]);
    assert(tmax[d] <= lsh[d]);
  }
}

GridPtrDesc::GridPtrDesc(const GHExt::LevelData &leveldata,
                         const MFPointer &mfp)
    : GridDesc(leveldata, mfp) {
  const amrex::Box &fbx = mfp.fabbox(); // allocated array
  cactus_offset = lbound(fbx);
}

GridPtrDesc1::GridPtrDesc1(const GHExt::LevelData::GroupData &groupdata,
                           const MFPointer &mfp)
    : GridDesc(groupdata.leveldata(), mfp) {
  const amrex::Box &fbx = mfp.fabbox(); // allocated array
  cactus_offset = lbound(fbx);
  for (int d = 0; d < dim; ++d) {
    assert(groupdata.nghostzones[d] >= 0);
    assert(groupdata.nghostzones[d] <= nghostzones[d]);
  }
  for (int d = 0; d < dim; ++d)
    gimin[d] = nghostzones[d] - groupdata.nghostzones[d];
  for (int d = 0; d < dim; ++d)
    gimax[d] = lsh[d] + (1 - groupdata.indextype[d]) -
               (nghostzones[d] - groupdata.nghostzones[d]);
  for (int d = 0; d < dim; ++d)
    gash[d] = ash[d] + (1 - groupdata.indextype[d]) -
              2 * (nghostzones[d] - groupdata.nghostzones[d]);
}

////////////////////////////////////////////////////////////////////////////////

// Ensure grid functions are valid
void error_if_invalid(const GHExt::LevelData::GroupData &groupdata, int vi,
                      int tl, const valid_t &required,
                      const function<string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VERROR("%s: Grid function \"%s\" is invalid on refinement level %d, "
                "time level %d; required %s, found %s",
                msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
                groupdata.level, tl, string(required).c_str(),
                string(groupdata.valid.at(tl).at(vi)).c_str());
}
void warn_if_invalid(const GHExt::LevelData::GroupData &groupdata, int vi,
                     int tl, const valid_t &required,
                     const function<string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VWARN(CCTK_WARN_ALERT,
               "%s: Grid function \"%s\" is invalid on refinement level %d, "
               "time level %d; required %s, found %s",
               msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
               groupdata.level, tl, string(required).c_str(),
               string(groupdata.valid.at(tl).at(vi)).c_str());
}

// Set grid functions to nan
void poison_invalid(const GHExt::LevelData::GroupData &groupdata, int vi,
                    int tl) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = groupdata.valid.at(tl).at(vi).get();
  if (valid.valid_all())
    return;

  const auto &leveldata = groupdata.leveldata();
  const auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling(
      {max_tile_size_x, max_tile_size_y, max_tile_size_z});
#pragma omp parallel
  for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {
    const GridPtrDesc1 grid(groupdata, mfi);
    const amrex::Array4<CCTK_REAL> &vars = groupdata.mfab.at(tl)->array(mfi);
    const GF3D1<CCTK_REAL> ptr_ = grid.gf3d(vars, vi);

    if (!valid.valid_any()) {
      grid.loop_idx(where_t::everywhere, groupdata.indextype,
                    groupdata.nghostzones,
                    [&](const Loop::PointDesc &p) { ptr_(p.I) = 0.0 / 0.0; });
    } else {
      if (!valid.valid_int)
        grid.loop_idx(where_t::interior, groupdata.indextype,
                      groupdata.nghostzones,
                      [&](const Loop::PointDesc &p) { ptr_(p.I) = 0.0 / 0.0; });
      if (!valid.valid_outer)
        grid.loop_idx(where_t::boundary, groupdata.indextype,
                      groupdata.nghostzones,
                      [&](const Loop::PointDesc &p) { ptr_(p.I) = 0.0 / 0.0; });
      if (!valid.valid_ghosts)
        grid.loop_idx(where_t::ghosts, groupdata.indextype,
                      groupdata.nghostzones,
                      [&](const Loop::PointDesc &p) { ptr_(p.I) = 0.0 / 0.0; });
    }
  }
}

// Ensure grid functions are not nan
// TODO: Parallelize this loop (and poison_invalid) further out, at
// the loop over time levels and grid variables
void check_valid(const GHExt::LevelData::GroupData &groupdata, int vi, int tl,
                 const function<string()> &msg) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = groupdata.valid.at(tl).at(vi).get();
  if (!valid.valid_any())
    return;

  size_t nan_count{0};
  array<int, 3> nan_imin, nan_imax;
  array<CCTK_REAL, 3> nan_xmin, nan_xmax;
  for (int d = 0; d < 3; ++d) {
    nan_imin[d] = numeric_limits<int>::max();
    nan_imax[d] = numeric_limits<int>::min();
    nan_xmin[d] = +1.0 / 0.0;
    nan_xmax[d] = -1.0 / 0.0;
  }
  const auto nan_update{
      [&](const GridDescBase &grid, const Loop::PointDesc &p) {
#pragma omp critical
        {
          ++nan_count;
          nan_imin[0] = min(nan_imin[0], grid.lbnd[0] + p.i);
          nan_imin[1] = min(nan_imin[1], grid.lbnd[1] + p.j);
          nan_imin[2] = min(nan_imin[2], grid.lbnd[2] + p.k);
          nan_imax[0] = max(nan_imax[0], grid.lbnd[0] + p.i);
          nan_imax[1] = max(nan_imax[1], grid.lbnd[1] + p.j);
          nan_imax[2] = max(nan_imax[2], grid.lbnd[2] + p.k);
          nan_xmin[0] = fmin(nan_xmin[0], p.x);
          nan_xmin[1] = fmin(nan_xmin[1], p.y);
          nan_xmin[2] = fmin(nan_xmin[2], p.z);
          nan_xmax[0] = fmax(nan_xmax[0], p.x);
          nan_xmax[1] = fmax(nan_xmax[1], p.y);
          nan_xmax[2] = fmax(nan_xmax[2], p.z);
        }
      }};
  const auto nan_check{[&](const GridDescBase &grid,
                           const GF3D1<const CCTK_REAL> &ptr_,
                           const Loop::PointDesc &p) {
    if (CCTK_BUILTIN_EXPECT(!CCTK_isfinite(ptr_(p.I)), false))
      nan_update(grid, p);
  }};
  const auto &leveldata = groupdata.leveldata();
  const auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling(
      {max_tile_size_x, max_tile_size_y, max_tile_size_z});
#pragma omp parallel
  for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {
    const GridPtrDesc1 grid(groupdata, mfi);
    const amrex::Array4<const CCTK_REAL> &vars =
        groupdata.mfab.at(tl)->array(mfi);
    const GF3D1<const CCTK_REAL> ptr_ = grid.gf3d(vars, vi);

    if (valid.valid_all()) {
      grid.loop_idx(
          where_t::everywhere, groupdata.indextype, groupdata.nghostzones,
          [&](const Loop::PointDesc &p) { nan_check(grid, ptr_, p); });
    } else {
      if (valid.valid_int)
        grid.loop_idx(
            where_t::interior, groupdata.indextype, groupdata.nghostzones,
            [&](const Loop::PointDesc &p) { nan_check(grid, ptr_, p); });
      if (valid.valid_outer)
        grid.loop_idx(
            where_t::boundary, groupdata.indextype, groupdata.nghostzones,
            [&](const Loop::PointDesc &p) { nan_check(grid, ptr_, p); });
      if (valid.valid_ghosts)
        grid.loop_idx(
            where_t::ghosts, groupdata.indextype, groupdata.nghostzones,
            [&](const Loop::PointDesc &p) { nan_check(grid, ptr_, p); });
    }
  }

  if (CCTK_BUILTIN_EXPECT(nan_count > 0, false)) {
#pragma omp critical
    {
      CCTK_VINFO(
          "%s: Grid function \"%s\" has %td nans on refinement level %d, time "
          "level %d, in box [%d,%d,%d]:[%d,%d,%d] (%g,%g,%g):(%g,%g,%g); "
          "expected valid %s",
          msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
          size_t(nan_count), leveldata.level, tl, nan_imin[0], nan_imin[1],
          nan_imin[2], nan_imax[0], nan_imax[1], nan_imax[2],
          double(nan_xmin[0]), double(nan_xmin[1]), double(nan_xmin[2]),
          double(nan_xmax[0]), double(nan_xmax[1]), double(nan_xmax[2]),
          string(groupdata.valid.at(tl).at(vi)).c_str());

      for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {
        const GridPtrDesc1 grid(groupdata, mfi);
        const amrex::Array4<const CCTK_REAL> &vars =
            groupdata.mfab.at(tl)->array(mfi);
        const GF3D1<const CCTK_REAL> ptr_ = grid.gf3d(vars, vi);

        if (valid.valid_int)
          grid.loop_idx(
              where_t::interior, groupdata.indextype, groupdata.nghostzones,
              [&](const Loop::PointDesc &p) {
                if (CCTK_BUILTIN_EXPECT(!CCTK_isfinite(ptr_(p.I)), false))
                  CCTK_VINFO("[%d,%d,%d] (%g,%g,%g) %g", p.i, p.j, p.k,
                             double(p.x), double(p.y), double(p.z),
                             double(ptr_(p.I)));
              });
        if (valid.valid_outer)
          grid.loop_idx(
              where_t::boundary, groupdata.indextype, groupdata.nghostzones,
              [&](const Loop::PointDesc &p) {
                if (CCTK_BUILTIN_EXPECT(!CCTK_isfinite(ptr_(p.I)), false))
                  CCTK_VINFO("[%d,%d,%d] (%g,%g,%g) %g", p.i, p.j, p.k,
                             double(p.x), double(p.y), double(p.z),
                             double(ptr_(p.I)));
              });
        if (valid.valid_ghosts)
          grid.loop_idx(
              where_t::ghosts, groupdata.indextype, groupdata.nghostzones,
              [&](const Loop::PointDesc &p) {
                if (CCTK_BUILTIN_EXPECT(!CCTK_isfinite(ptr_(p.I)), false))
                  CCTK_VINFO("[%d,%d,%d] (%g,%g,%g) %g", p.i, p.j, p.k,
                             double(p.x), double(p.y), double(p.z),
                             double(ptr_(p.I)));
              });
      }

      CCTK_VERROR(
          "%s: Grid function \"%s\" has nans on refinement level %d, time "
          "level %d; expected valid %s",
          msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
          leveldata.level, tl, string(groupdata.valid.at(tl).at(vi)).c_str());
    }
  }
}

// Ensure grid functions are valid
void error_if_invalid(const GHExt::GlobalData::ScalarGroupData &groupdata,
                      int vi, int tl, const valid_t &required,
                      const function<string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VERROR("%s: Grid function \"%s\" is invalid on time level %d; "
                "required %s, found %s",
                msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
                tl, string(required).c_str(),
                string(groupdata.valid.at(tl).at(vi)).c_str());
}
void warn_if_invalid(const GHExt::GlobalData::ScalarGroupData &groupdata,
                     int vi, int tl, const valid_t &required,
                     const function<string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VWARN(CCTK_WARN_ALERT,
               "%s: Grid function \"%s\" is invalid on time level %d; "
               "required %s, found %s",
               msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
               tl, string(required).c_str(),
               string(groupdata.valid.at(tl).at(vi)).c_str());
}

// Set grid scalars to nan
void poison_invalid(const GHExt::GlobalData::ScalarGroupData &scalargroupdata,
                    int vi, int tl) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = scalargroupdata.valid.at(tl).at(vi).get();
  if (valid.valid_all())
    return;

  // scalars have no boundary so we expect them to alway be valid
  assert(valid.valid_outer && valid.valid_ghosts);

  if (!valid.valid_int) {
    CCTK_REAL *restrict const ptr =
        const_cast<CCTK_REAL *>(&scalargroupdata.data.at(tl).at(vi));
    *ptr = 0.0 / 0.0;
  }
}

// Ensure grid scalars are not nan
void check_valid(const GHExt::GlobalData::ScalarGroupData &scalargroupdata,
                 int vi, int tl, const function<string()> &msg) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = scalargroupdata.valid.at(tl).at(vi).get();
  if (!valid.valid_any())
    return;

  // scalars have no boundary so we expect them to alway be valid
  assert(valid.valid_outer && valid.valid_ghosts);

  atomic<size_t> nan_count{0};
  if (valid.valid_int) {
    const CCTK_REAL *restrict const ptr = &scalargroupdata.data.at(tl).at(vi);
    if (CCTK_BUILTIN_EXPECT(!CCTK_isfinite(*ptr), false)) {
      ++nan_count;
    }
  }

  if (CCTK_BUILTIN_EXPECT(nan_count > 0, false))
    CCTK_VERROR("%s: Grid Scalar \"%s\" has %td nans on time level %d; "
                "expected valid %s",
                msg().c_str(),
                CCTK_FullVarName(scalargroupdata.firstvarindex + vi),
                size_t(nan_count), tl,
                string(scalargroupdata.valid.at(tl).at(vi)).c_str());
}

struct tiletag_t {
  int level;
  amrex::Box tilebox;
  int gi, vi, tl;
  tiletag_t() = delete;

  friend bool operator==(const tiletag_t &x, const tiletag_t &y) {
    return make_tuple(x.level, x.tilebox, x.gi, x.vi, x.tl) ==
           make_tuple(y.level, y.tilebox, y.gi, y.vi, y.tl);
  }
  friend bool operator<(const tiletag_t &x, const tiletag_t &y) {
    return make_tuple(x.level, x.tilebox, x.gi, x.vi, x.tl) <
           make_tuple(y.level, y.tilebox, y.gi, y.vi, y.tl);
  }

  friend ostream &operator<<(ostream &os, const tiletag_t &x) {
    return os << "tiletag_t{"
              << "level:" << x.level << ","
              << "tilebox:" << x.tilebox << ","
              << "gi:" << x.gi << ","
              << "vi:" << x.vi << ","
              << "tl:" << x.tl << "}";
  }
  operator string() const {
    ostringstream buf;
    buf << *this;
    return buf.str();
  }
};

struct checksum_t {
  valid_t where;
  uLong crc;
  checksum_t() = default;
  inline checksum_t(const valid_t &where)
      : where(where), crc(crc32(0, nullptr, 0)) {}
  template <typename T> inline void add(const T &x) {
    crc = crc32(crc, static_cast<const Bytef *>(static_cast<const void *>(&x)),
                sizeof x);
  }

  friend bool operator==(const checksum_t &x, const checksum_t &y) {
    return x.where == y.where && x.crc == y.crc;
  }
  friend bool operator!=(const checksum_t &x, const checksum_t &y) {
    return !(x == y);
  }

  friend ostream &operator<<(ostream &os, const checksum_t &x) {
    return os << "checksum_t{where:" << x.where << ",crc:0x" << hex
              << setfill('0') << setw(8) << x.crc << "}";
  }
  operator string() const {
    ostringstream buf;
    buf << *this;
    return buf.str();
  }
};

typedef map<tiletag_t, checksum_t> checksums_t;

checksums_t
calculate_checksums(const vector<vector<vector<valid_t> > > &will_write) {
  DECLARE_CCTK_PARAMETERS;

  checksums_t checksums;

  if (!poison_undefined_values)
    return checksums;

  assert(active_levels);
  active_levels->loop([&](auto &restrict leveldata) {
    auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling(
        {max_tile_size_x, max_tile_size_y, max_tile_size_z});
#pragma omp parallel
    for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {

      for (const auto &groupdataptr : leveldata.groupdata) {
        if (groupdataptr == nullptr)
          continue;

        auto &restrict groupdata = *groupdataptr;
        const GridPtrDesc1 grid(groupdata, mfi);

        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          for (int tl = 0; tl < int(groupdata.valid.size()); ++tl) {
            const tiletag_t tiletag{leveldata.level, mfi.tilebox(),
                                    groupdata.groupindex, vi, tl};

            const auto &valid = groupdata.valid.at(tl).at(vi).get();
            // No information given for this timelevel; assume not written
            if (tl >= int(will_write.at(groupdata.groupindex).at(vi).size()))
              continue;
            const auto &wr = will_write.at(groupdata.groupindex).at(vi).at(tl);
            valid_t to_check = valid & ~wr;

            // Check only those variables which are valid, and where
            // some part (but not everything) is written
            if (!(wr.valid_any() && to_check.valid_any()))
              continue;

            const amrex::Array4<const CCTK_REAL> &vars =
                groupdata.mfab.at(tl)->array(mfi);
            const GF3D1<const CCTK_REAL> var_ = grid.gf3d(vars, vi);

            checksum_t checksum(to_check);
            checksum.add(tiletag);
            const auto add_point{
                [&](const Loop::PointDesc &p) { checksum.add(var_(p.I)); }};

            if (to_check.valid_int)
              grid.loop_idx(where_t::interior, groupdata.indextype,
                            groupdata.nghostzones, add_point);

            if (to_check.valid_outer)
              grid.loop_idx(where_t::boundary, groupdata.indextype,
                            groupdata.nghostzones, add_point);

            if (to_check.valid_ghosts)
              grid.loop_idx(where_t::ghosts, groupdata.indextype,
                            groupdata.nghostzones, add_point);

#pragma omp critical
            checksums[tiletag] = checksum;
          }
        }
      }
    }
  });
  return checksums;
}

void check_checksums(const checksums_t checksums) {
  DECLARE_CCTK_PARAMETERS;

  if (!poison_undefined_values)
    return;
  if (checksums.empty())
    return;

  assert(active_levels);
  active_levels->loop([&](auto &restrict leveldata) {
    auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling(
        {max_tile_size_x, max_tile_size_y, max_tile_size_z});
#pragma omp parallel
    for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {

      for (const auto &groupdataptr : leveldata.groupdata) {
        if (groupdataptr == nullptr)
          continue;

        auto &restrict groupdata = *groupdataptr;
        const GridPtrDesc1 grid(groupdata, mfi);

        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          for (int tl = 0; tl < int(groupdata.valid.size()); ++tl) {
            const tiletag_t tiletag{leveldata.level, mfi.tilebox(),
                                    groupdata.groupindex, vi, tl};

            if (!checksums.count(tiletag))
              continue;

            const auto &old_checksum = checksums.at(tiletag);
            const auto &did_check = old_checksum.where;
            assert(did_check.valid_any());

            const amrex::Array4<const CCTK_REAL> &vars =
                groupdata.mfab.at(tl)->array(mfi);
            const GF3D1<const CCTK_REAL> var_ = grid.gf3d(vars, vi);

            checksum_t checksum(did_check);
            checksum.add(tiletag);
            const auto add_point{
                [&](const Loop::PointDesc &p) { checksum.add(var_(p.I)); }};

            if (did_check.valid_int)
              grid.loop_idx(where_t::interior, groupdata.indextype,
                            groupdata.nghostzones, add_point);

            if (did_check.valid_outer)
              grid.loop_idx(where_t::boundary, groupdata.indextype,
                            groupdata.nghostzones, add_point);

            if (did_check.valid_ghosts)
              grid.loop_idx(where_t::ghosts, groupdata.indextype,
                            groupdata.nghostzones, add_point);

            if (checksum != old_checksum)
              CCTK_VERROR(
                  "Checksum mismatch: variable %s, tile %s, "
                  "int:%d,outer:%d,ghosts:%d, old checksum %s, new checksum %s",
                  CCTK_FullVarName(groupdata.firstvarindex + tiletag.vi),
                  string(tiletag).c_str(), int(did_check.valid_int),
                  int(did_check.valid_outer), int(did_check.valid_ghosts),
                  string(old_checksum).c_str(), string(checksum).c_str());
          }
        }
      }
    }
  });
}

////////////////////////////////////////////////////////////////////////////////

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

enum class mode_t { unknown, local, level, global, meta };

mode_t current_mode(const cGH *restrict cctkGH) {
  if (cctkGH->cctk_lsh[0] != undefined)
    return mode_t::local;
  else if (cctkGH->cctk_gsh[0] != undefined)
    return mode_t::level;
  else if (cctkGH->cctk_nghostzones[0] != undefined)
    return mode_t::global;
  else
    return mode_t::meta;
}

bool in_local_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::local;
}

bool in_level_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::level;
}

bool in_global_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::global;
}

bool in_meta_mode(const cGH *restrict cctkGH) {
  return current_mode(cctkGH) == mode_t::meta;
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
  const amrex::Geometry &geom = ghext->amrcore->Geom(0);
  const CCTK_REAL *restrict const x0 = geom.ProbLo();
  const CCTK_REAL *restrict const dx = geom.CellSize();

  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_origin_space[d] = NAN;
    cctkGH->cctk_delta_space[d] = NAN;
  }

  // Initialize time stepping
  // CCTK_REAL mindx = 1.0 / 0.0;
  // const int numlevels = ghext->amrcore->finestLevel() + 1;
  // for (int level = 0; level < numlevels; ++level) {
  //   const amrex::Geometry &geom = ghext->amrcore->Geom(level);
  //   const CCTK_REAL *restrict dx = geom.CellSize();
  //   for (int d = 0; d < dim; ++d)
  //     mindx = fmin(mindx, dx[d]);
  // }
  CCTK_REAL mindx = 1.0 / 0.0;
  for (int d = 0; d < dim; ++d)
    mindx = fmin(mindx, dx[d]);
  mindx /= 1 << (max_num_levels - 1);
  cctkGH->cctk_time = 0.0;
  cctkGH->cctk_delta_time = dtfac * mindx;
  // init into meta mode
  cctkGH->cctk_nghostzones[0] = undefined;
  cctkGH->cctk_lsh[0] = undefined;
  cctkGH->cctk_gsh[0] = undefined;
  assert(in_meta_mode(cctkGH));
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
  assert(in_meta_mode(cctkGH));

  // The number of ghostzones in each direction
  // TODO: Get this from mfab (mfab.fb_ghosts)
  cctkGH->cctk_nghostzones[0] = ghost_size >= 0 ? ghost_size : ghost_size_x;
  cctkGH->cctk_nghostzones[1] = ghost_size >= 0 ? ghost_size : ghost_size_y;
  cctkGH->cctk_nghostzones[2] = ghost_size >= 0 ? ghost_size : ghost_size_z;

  // Grid scalar pointers
  {
    auto &restrict globaldata = ghext->globaldata;
    const int num_groups = CCTK_NumGroups();
    for (int gi = 0; gi < num_groups; ++gi) {
      cGroup group;
      int ierr = CCTK_GroupData(gi, &group);
      assert(!ierr);

      if (group.grouptype != CCTK_SCALAR)
        continue;

      auto &restrict scalargroupdata = *globaldata.scalargroupdata.at(gi);
      for (int tl = 0; tl < int(scalargroupdata.data.size()); ++tl) {
        const auto &restrict vars = scalargroupdata.data.at(tl);
        for (int vi = 0; vi < scalargroupdata.numvars; ++vi) {
          cctkGH->data[scalargroupdata.firstvarindex + vi][tl] =
              const_cast<CCTK_REAL *>(&vars.at(vi));
        }
      }
    }
  }
}
void leave_global_mode(cGH *restrict cctkGH) {
  assert(in_global_mode(cctkGH));

  // Grid scalar pointers
  {
    auto &restrict globaldata = ghext->globaldata;
    const int num_groups = CCTK_NumGroups();
    for (int gi = 0; gi < num_groups; ++gi) {
      cGroup group;
      int ierr = CCTK_GroupData(gi, &group);
      assert(!ierr);

      if (group.grouptype != CCTK_SCALAR)
        continue;

      auto &restrict scalargroupdata = *globaldata.scalargroupdata.at(gi);
      for (int tl = 0; tl < int(scalargroupdata.data.size()); ++tl) {
        for (int vi = 0; vi < scalargroupdata.numvars; ++vi) {
          cctkGH->data[scalargroupdata.firstvarindex + vi][tl] = nullptr;
        }
      }
    }
  }

  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_nghostzones[d] = undefined;
}

// Set cctkGH entries for level mode
void enter_level_mode(cGH *restrict cctkGH,
                      const GHExt::LevelData &restrict leveldata) {
  DECLARE_CCTK_PARAMETERS;
  assert(in_global_mode(cctkGH));

  // Global shape
  const amrex::Box &domain = ghext->amrcore->Geom(leveldata.level).Domain();
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_gsh[d] = domain[orient(d, 1)] - domain[orient(d, 0)] + 1 +
                          2 * cctkGH->cctk_nghostzones[d];

  const amrex::Geometry &geom = ghext->amrcore->Geom(0);
  const CCTK_REAL *restrict const global_x0 = geom.ProbLo();
  const CCTK_REAL *restrict const global_dx = geom.CellSize();
  for (int d = 0; d < dim; ++d) {
    // The refinement factor over the top level (coarsest) grid
    const int levfac = 1 << leveldata.level;
    cctkGH->cctk_levfac[d] = levfac;
    // Offset between this level's and the coarsest level's origin as multiple
    // of the grid spacing
    const int levoff = 1 - 2 * cctkGH->cctk_nghostzones[d];
    const int levoffdenom = 2;
    cctkGH->cctk_levoff[d] = levoff;
    cctkGH->cctk_levoffdenom[d] = levoffdenom;
    // Coordinates
    const CCTK_REAL origin_space =
        global_x0[d] +
        (1 - cctkGH->cctk_nghostzones[d] + CCTK_REAL(1) / 2) * global_dx[d];
    const CCTK_REAL delta_space = global_dx[d];
    cctkGH->cctk_delta_space[d] = delta_space / levfac;
    cctkGH->cctk_origin_space[d] =
        origin_space + cctkGH->cctk_delta_space[d] * levoff / levoffdenom;
  }
}
void leave_level_mode(cGH *restrict cctkGH,
                      const GHExt::LevelData &restrict leveldata) {
  assert(in_level_mode(cctkGH));
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_gsh[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_levfac[d] = undefined;
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_levoff[d] = undefined;
    cctkGH->cctk_levoffdenom[d] = 0;
  }
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_origin_space[d] = NAN;
    cctkGH->cctk_delta_space[d] = NAN;
  }
}

// Set cctkGH entries for local mode
// TODO: Have separate cctkGH for each level, each local box, and each tile
void enter_local_mode(cGH *restrict cctkGH, TileBox &restrict tilebox,
                      const GHExt::LevelData &restrict leveldata,
                      const MFPointer &mfp) {
  assert(in_level_mode(cctkGH));
  const GridPtrDesc grid(leveldata, mfp);

  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lsh[d] = grid.lsh[d];
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ash[d] = grid.ash[d];
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_lbnd[d] = grid.lbnd[d];
    cctkGH->cctk_ubnd[d] = grid.ubnd[d];
  }
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      cctkGH->cctk_bbox[2 * d + f] = grid.bbox[2 * d + f];

  for (int d = 0; d < dim; ++d) {
    tilebox.tile_min[d] = grid.tmin[d];
    tilebox.tile_max[d] = grid.tmax[d];
  }

  // Grid function pointers
  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata = *leveldata.groupdata.at(gi);
    const GridPtrDesc1 grid1(groupdata, mfp);
    for (int tl = 0; tl < int(groupdata.mfab.size()); ++tl) {
      const amrex::Array4<CCTK_REAL> vars =
          groupdata.mfab.at(tl)->array(mfp.index());
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        cctkGH->data[groupdata.firstvarindex + vi][tl] = grid1.ptr(vars, vi);
    }
  }

  // Check constraints
  for (int d = 0; d < dim; ++d) {
    // Domain size
    assert(cctkGH->cctk_gsh[d] >= 0);

    // Local size
    assert(cctkGH->cctk_lbnd[d] >= 0);
    assert(cctkGH->cctk_lsh[d] >= 0);
    assert(cctkGH->cctk_lbnd[d] + cctkGH->cctk_lsh[d] <= cctkGH->cctk_gsh[d]);
    assert(cctkGH->cctk_ubnd[d] ==
           cctkGH->cctk_lbnd[d] + cctkGH->cctk_lsh[d] - 1);

    // Internal representation
    assert(cctkGH->cctk_ash[d] >= 0);
    assert(cctkGH->cctk_ash[d] >= cctkGH->cctk_lsh[d]);

    // Ghost zones
    assert(cctkGH->cctk_nghostzones[d] >= 0);
    assert(2 * cctkGH->cctk_nghostzones[d] <= cctkGH->cctk_lsh[d]);

    // Tiles
    assert(tilebox.tile_min[d] >= 0);
    assert(tilebox.tile_max[d] >= tilebox.tile_min[d]);
    assert(tilebox.tile_max[d] <= cctkGH->cctk_lsh[d]);
  }
}
void leave_local_mode(cGH *restrict cctkGH, TileBox &restrict tilebox,
                      const GHExt::LevelData &restrict leveldata,
                      const MFPointer &mfp) {
  assert(in_local_mode(cctkGH));
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_lsh[d] = undefined;
  for (int d = 0; d < dim; ++d)
    cctkGH->cctk_ash[d] = undefined;
  for (int d = 0; d < dim; ++d) {
    cctkGH->cctk_lbnd[d] = undefined;
    cctkGH->cctk_ubnd[d] = undefined;
  }
  for (int d = 0; d < dim; ++d)
    for (int f = 0; f < 2; ++f)
      cctkGH->cctk_bbox[2 * d + f] = undefined;
  for (int d = 0; d < dim; ++d) {
    tilebox.tile_min[d] = undefined;
    tilebox.tile_max[d] = undefined;
  }
  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata = *leveldata.groupdata.at(gi);
    for (int tl = 0; tl < int(groupdata.mfab.size()); ++tl) {
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        cctkGH->data[groupdata.firstvarindex + vi][tl] = nullptr;
    }
  }
}
extern "C" void CarpetX_GetTileExtent(const void *restrict cctkGH_,
                                      CCTK_INT *restrict tile_min,
                                      CCTK_INT *restrict tile_max) {
  const cGH *restrict cctkGH = static_cast<const cGH *>(cctkGH_);
  // Check whether we are in local mode
  assert(cctkGH->cctk_bbox[0] != undefined);
  int thread_num = omp_get_thread_num();
  if (omp_in_parallel()) {
    const cGH *restrict threadGH = &thread_local_info.at(thread_num)->cctkGH;
    // Check whether this is the correct cGH structure
    assert(cctkGH == threadGH);
  }

  const TileBox &restrict tilebox = thread_local_info.at(thread_num)->tilebox;
  for (int d = 0; d < dim; ++d) {
    tile_min[d] = tilebox.tile_min[d];
    tile_max[d] = tilebox.tile_max[d];
  }
}

extern "C" void CarpetX_GetLoopBoxAll(const void *restrict const cctkGH_,
                                      const CCTK_INT size,
                                      CCTK_INT *restrict const loop_min,
                                      CCTK_INT *restrict const loop_max) {
  const cGH *restrict const cctkGH = static_cast<const cGH *>(cctkGH_);
  assert(size == dim);
  const GridDescBase grid(cctkGH);
  vect<int, dim> imin, imax;
  // Assumd cell centred box
  grid.box_all<1, 1, 1>(grid.nghostzones, imin, imax);
  for (int d = 0; d < dim; ++d)
    loop_min[d] = imin[d];
  for (int d = 0; d < dim; ++d)
    loop_max[d] = imax[d];
}

extern "C" void CarpetX_GetLoopBoxInt(const void *restrict const cctkGH_,
                                      const CCTK_INT size,
                                      CCTK_INT *restrict const loop_min,
                                      CCTK_INT *restrict const loop_max) {
  const cGH *restrict const cctkGH = static_cast<const cGH *>(cctkGH_);
  assert(size == dim);
  const GridDescBase grid(cctkGH);
  vect<int, dim> imin, imax;
  // Assumd cell centred box
  grid.box_int<1, 1, 1>(grid.nghostzones, imin, imax);
  for (int d = 0; d < dim; ++d)
    loop_min[d] = imin[d];
  for (int d = 0; d < dim; ++d)
    loop_max[d] = imax[d];
}

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

enum class rdwr_t { read, write, invalid };
ostream &operator<<(ostream &os, const rdwr_t rdwr) {
  switch (rdwr) {
  case rdwr_t::read:
    return os << "read";
  case rdwr_t::write:
    return os << "write";
  case rdwr_t::invalid:
    return os << "invalid";
  default:
    assert(0);
  }
}

struct clause_t {
  int gi, vi, tl;
  valid_t valid;

  friend bool operator==(const clause_t &x, const clause_t &y) {
    return make_tuple(x.gi, x.vi, x.tl, x.valid) ==
           make_tuple(y.gi, y.vi, y.tl, y.valid);
  }
  friend bool operator<(const clause_t &x, const clause_t &y) {
    return make_tuple(x.gi, x.vi, x.tl, x.valid) <
           make_tuple(y.gi, y.vi, y.tl, y.valid);
  }

  friend ostream &operator<<(ostream &os, const clause_t &cl) {
    return os << "clause_t{gi:" << cl.gi << ",vi:" << cl.vi << ",tl:" << cl.tl
              << ",valid:" << cl.valid << "}";
  }
};

vector<clause_t> decode_clauses(const cFunctionData *restrict attribute,
                                const rdwr_t rdwr) {
  vector<clause_t> result;
  result.reserve(attribute->n_RDWR);
  for (int n = 0; n < attribute->n_RDWR; ++n) {
    const RDWR_entry &restrict RDWR = attribute->RDWR[n];
    int gi = CCTK_GroupIndexFromVarI(RDWR.varindex);
    assert(gi >= 0);
    int vi = RDWR.varindex - CCTK_FirstVarIndexI(gi);
    assert(vi >= 0 && vi < CCTK_NumVarsInGroupI(gi));
    int tl = RDWR.timelevel;
    int where;
    switch (rdwr) {
    case rdwr_t::read:
      where = RDWR.where_rd;
      break;
    case rdwr_t::write:
      where = RDWR.where_wr;
      break;
    case rdwr_t::invalid:
      where = RDWR.where_inv;
      break;
    default:
      assert(0);
    }
    valid_t valid;
    valid.valid_int = where & CCTK_VALID_INTERIOR;
    valid.valid_outer = where & CCTK_VALID_BOUNDARY;
    valid.valid_ghosts = where & CCTK_VALID_GHOSTS;
    result.push_back({gi, vi, tl, valid});
  }
  return result;
}

// Schedule initialisation
int Initialise(tFleshConfig *config) {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("Initialise");
  Interval interval(timer);

  cGH *restrict const cctkGH = CCTK_SetupGH(config, 0);
  CCTKi_AddGH(config, 0, cctkGH);

  // Initialise iteration and time
  cctkGH->cctk_iteration = 0;
  cctkGH->cctk_time = *static_cast<const CCTK_REAL *>(
      CCTK_ParameterGet("cctk_initial_time", "Cactus", nullptr));

  // Initialise schedule
  CCTKi_ScheduleGHInit(cctkGH);

  // Check presync mode
  if (!CCTK_EQUALS(presync_mode, "mixed-error"))
    CCTK_ERROR(
        "CarpetX currently requires Cactus::presync_mode = \"mixed-error\"");

  // Initialise all grid extensions
  CCTKi_InitGHExtensions(cctkGH);

  // Set up cctkGH
  setup_cctkGH(cctkGH);
  enter_global_mode(cctkGH);

  int max_threads = omp_get_max_threads();
  thread_local_info.resize(max_threads);
  for (int n = 0; n < max_threads; ++n) {
    thread_local_info.at(n) = make_unique<thread_local_info_t>();
    cGH *restrict threadGH = &thread_local_info.at(n)->cctkGH;
    clone_cctkGH(threadGH, cctkGH);
    setup_cctkGH(threadGH);
    enter_global_mode(threadGH);
    thread_local_info.at(n)->mfpointer = nullptr;
  }
  swap(saved_thread_local_info, thread_local_info);
  assert(thread_local_info.empty());

  assert(ghext->leveldata.empty());
  assert(!active_levels);
  active_levels = make_optional<active_levels_t>(0, 0);

  CCTK_Traverse(cctkGH, "CCTK_WRAGH");
  CCTK_Traverse(cctkGH, "CCTK_PARAMCHECK");
  CCTKi_FinaliseParamWarn();

  active_levels = optional<active_levels_t>();

  if (config->recovered) {
    // Recover
#pragma omp critical
    CCTK_VINFO("Recovering from checkpoint...");

    const char *recovery_mode = *static_cast<const char *const *>(
        CCTK_ParameterGet("recovery_mode", "Cactus", nullptr));

    CCTK_Traverse(cctkGH, "CCTK_BASEGRID");

    if (!CCTK_Equals(recovery_mode, "strict")) {
      // Set up initial conditions
      CCTK_Traverse(cctkGH, "CCTK_INITIAL");
      CCTK_Traverse(cctkGH, "CCTK_POSTINITIAL");
      CCTK_Traverse(cctkGH, "CCTK_POSTPOSTINITIAL");
    }

    // Recover
    CCTK_Traverse(cctkGH, "CCTK_RECOVER_VARIABLES");
    CCTK_Traverse(cctkGH, "CCTK_POST_RECOVER_VARIABLES");

  } else {
    // Set up initial conditions
#pragma omp critical
    CCTK_VINFO("Setting up initial conditions...");

    // Create coarse grid
    {
      static Timer timer("InitialiseRegrid");
      Interval interval(timer);
      const CCTK_REAL time = 0.0; // dummy time
      assert(saved_cctkGH == nullptr);
      saved_cctkGH = cctkGH;
      ghext->amrcore->MakeNewGrids(time);
      saved_cctkGH = nullptr;
    }

    // Output domain information
    if (CCTK_MyProc(nullptr) == 0) {
      enter_level_mode(cctkGH, ghext->leveldata.at(0));
      const int *restrict const gsh = cctkGH->cctk_gsh;
      const int *restrict const nghostzones = cctkGH->cctk_nghostzones;
      CCTK_REAL x0[dim], x1[dim], dx[dim];
      for (int d = 0; d < dim; ++d) {
        dx[d] = cctkGH->cctk_delta_space[d];
        x0[d] = cctkGH->cctk_origin_space[d] - dx[d] / 2;
        x1[d] = x0[d] + (gsh[d] - 2 * nghostzones[d]) * dx[d];
      }
#pragma omp critical
      {
        CCTK_VINFO("Grid extent:");
        CCTK_VINFO("  gsh=[%d,%d,%d]", gsh[0], gsh[1], gsh[2]);
        CCTK_VINFO("Domain extent:");
        CCTK_VINFO("  xmin=[%.17g,%.17g,%.17g]", double(x0[0]), double(x0[1]),
                   double(x0[2]));
        CCTK_VINFO("  xmax=[%.17g,%.17g,%.17g]", double(x1[0]), double(x1[1]),
                   double(x1[2]));
        CCTK_VINFO("  base dx=[%.17g,%.17g,%.17g]", double(dx[0]),
                   double(dx[1]), double(dx[2]));
        CCTK_VINFO("Time stepping:");
        CCTK_VINFO("  t0=%.17g", double(cctkGH->cctk_time));
        CCTK_VINFO("  dt=%.17g", double(cctkGH->cctk_delta_time));
      }
      leave_level_mode(cctkGH, ghext->leveldata.at(0));
    }

    for (;;) {
      const int level = ghext->amrcore->finestLevel();
#pragma omp critical
      CCTK_VINFO("Initializing level %d...", level);

      assert(!active_levels);
      active_levels = make_optional<active_levels_t>(0, level + 1);

      InputGH(cctkGH);
      CCTK_Traverse(cctkGH, "CCTK_INITIAL");
      CCTK_Traverse(cctkGH, "CCTK_POSTINITIAL");
      CCTK_Traverse(cctkGH, "CCTK_POSTPOSTINITIAL");

      active_levels = optional<active_levels_t>();

      if (level >= ghext->amrcore->maxLevel())
        break;

#pragma omp critical
      CCTK_VINFO("Regridding...");
      const int old_numlevels = ghext->amrcore->finestLevel() + 1;
      {
        static Timer timer("InitialiseRegrid");
        Interval interval(timer);
        assert(saved_cctkGH == nullptr);
        saved_cctkGH = cctkGH;
        const CCTK_REAL time = 0.0; // dummy time
        ghext->amrcore->regrid(0, time);
        saved_cctkGH = nullptr;
      }
      const int new_numlevels = ghext->amrcore->finestLevel() + 1;
      const int max_numlevels = ghext->amrcore->maxLevel() + 1;
      assert(new_numlevels >= 0 && new_numlevels <= max_numlevels);
      assert(new_numlevels == old_numlevels ||
             new_numlevels == old_numlevels + 1);
#pragma omp critical
      {
        const double pts0 = ghext->leveldata.at(0).fab->boxArray().d_numPts();
        for (const auto &leveldata : ghext->leveldata) {
          const int sz = leveldata.fab->size();
          const double pts = leveldata.fab->boxArray().d_numPts();
          if (leveldata.level == 0) {
            CCTK_VINFO("  level %d: %d boxes, %.0f cells (%.4g%%)",
                       leveldata.level, sz, pts,
                       100 * pts / (pow(2.0, dim * leveldata.level) * pts0));
          } else {
            const double ptsc = ghext->leveldata.at(leveldata.level - 1)
                                    .fab->boxArray()
                                    .d_numPts();
            CCTK_VINFO("  level %d: %d boxes, %.0f cells (%.4g%%, %.0f%%)",
                       leveldata.level, sz, pts,
                       100 * pts / (pow(2.0, dim * leveldata.level) * pts0),
                       100 * pts / (pow(2.0, dim) * ptsc));
          }
        }
      }

      // Did we create a new level?
      const bool did_create_new_level = new_numlevels > old_numlevels;
      if (!did_create_new_level)
        break;
    }
  }
#pragma omp critical
  CCTK_VINFO("Initialized %d levels", int(ghext->leveldata.size()));

  assert(!active_levels);
  active_levels = make_optional<active_levels_t>(0, ghext->leveldata.size());

  if (!restrict_during_sync) {
    // Restrict
    assert(active_levels);
    active_levels->loop_reverse([&](const auto &leveldata) {
      if (leveldata.level != int(ghext->leveldata.size()) - 1)
        Restrict(leveldata.level);
    });
    CCTK_Traverse(cctkGH, "CCTK_POSTRESTRICT");
  }

  // Checkpoint, analysis, output
  CCTK_Traverse(cctkGH, "CCTK_POSTSTEP");
  CCTK_Traverse(cctkGH, "CCTK_CPINITIAL");
  CCTK_Traverse(cctkGH, "CCTK_ANALYSIS");
  CCTK_OutputGH(cctkGH);

  active_levels = optional<active_levels_t>();

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

void InvalidateTimelevels(cGH *restrict const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype == CCTK_GF) {
      assert(active_levels);
      active_levels->loop([&](const auto &restrict leveldata) {
        auto &restrict groupdata = *leveldata.groupdata.at(gi);
        if (!groupdata.do_checkpoint) {
          // Invalidate all time levels
          const int ntls = groupdata.mfab.size();
          for (int tl = 0; tl < ntls; ++tl) {
            for (int vi = 0; vi < groupdata.numvars; ++vi) {
              groupdata.valid.at(tl).at(vi).set(valid_t(), [] {
                return "InvalidateTimelevels (invalidate all "
                       "non-checkpointed variables)";
              });
              poison_invalid(groupdata, vi, tl);
            }
          }
        }
      });
    }

    // Invalidate scalars
    if (group.grouptype == CCTK_SCALAR) {

      auto &restrict globaldata = ghext->globaldata;
      auto &restrict scalargroupdata = *globaldata.scalargroupdata.at(gi);
      if (!scalargroupdata.do_checkpoint) {
        // Invalidate all time levels
        const int ntls = scalargroupdata.data.size();
        for (int tl = 0; tl < ntls; ++tl) {
          for (int vi = 0; vi < scalargroupdata.numvars; ++vi) {
            // TODO: handle this more nicely
            scalargroupdata.valid.at(tl).at(vi).set_int(false, [] {
              return "InvalidateTimelevels (invalidate all non-checkpointed "
                     "variables)";
            });
            poison_invalid(scalargroupdata, vi, tl);
          }
        }
      }
    }

  } // for gi
}

void CycleTimelevels(cGH *restrict const cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  cctkGH->cctk_iteration += 1;
  cctkGH->cctk_time += cctkGH->cctk_delta_time;

  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype == CCTK_GF) {

      assert(active_levels);
      active_levels->loop([&](auto &restrict leveldata) {
        auto &restrict groupdata = *leveldata.groupdata.at(gi);
        const int ntls = groupdata.mfab.size();
        // Rotate time levels and invalidate current time level
        if (ntls > 1) {
          rotate(groupdata.mfab.begin(), groupdata.mfab.end() - 1,
                 groupdata.mfab.end());
          rotate(groupdata.valid.begin(), groupdata.valid.end() - 1,
                 groupdata.valid.end());
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            groupdata.valid.at(0).at(vi).set(valid_t(), [] {
              return "CycletimeLevels (invalidate current time level)";
            });
            poison_invalid(groupdata, vi, 0);
          }
        }
        for (int tl = 0; tl < ntls; ++tl)
          for (int vi = 0; vi < groupdata.numvars; ++vi)
            check_valid(groupdata, vi, tl, [&]() { return "CycleTimelevels"; });
      });
    }

    // cycle scalars
    if (group.grouptype == CCTK_SCALAR) {

      auto &restrict globaldata = ghext->globaldata;
      auto &restrict scalargroupdata = *globaldata.scalargroupdata.at(gi);
      const int ntls = scalargroupdata.data.size();
      // Rotate time levels and invalidate current time level
      if (ntls > 1) {
        rotate(scalargroupdata.data.begin(), scalargroupdata.data.end() - 1,
               scalargroupdata.data.end());
        rotate(scalargroupdata.valid.begin(), scalargroupdata.valid.end() - 1,
               scalargroupdata.valid.end());
        for (int vi = 0; vi < scalargroupdata.numvars; ++vi) {
          scalargroupdata.valid.at(0).at(vi).set_int(false, [] {
            return "CycletimeLevels (invalidate current time level)";
          });
          poison_invalid(scalargroupdata, vi, 0);
        }
      }
      for (int tl = 0; tl < ntls; ++tl)
        for (int vi = 0; vi < scalargroupdata.numvars; ++vi)
          check_valid(scalargroupdata, vi, tl,
                      [&]() { return "CycleTimelevels"; });
    }

  } // for gi
}

// Schedule evolution
int Evolve(tFleshConfig *config) {
  DECLARE_CCTK_PARAMETERS;

  static Timer timer("Evolve");
  Interval interval(timer);

  assert(config);
  cGH *restrict const cctkGH = config->GH[0];
  assert(cctkGH);

#pragma omp critical
  CCTK_VINFO("Starting evolution...");

  while (!EvolutionIsDone(cctkGH)) {

    assert(!active_levels);

    // TODO: Move regridding into a function
    if (regrid_every > 0 && cctkGH->cctk_iteration % regrid_every == 0 &&
        ghext->amrcore->maxLevel() > 0) {
#pragma omp critical
      CCTK_VINFO("Regridding...");
      const int old_numlevels = ghext->amrcore->finestLevel() + 1;
      {
        static Timer timer("EvolveRegrid");
        Interval interval(timer);
        assert(saved_cctkGH == nullptr);
        saved_cctkGH = cctkGH;
        CCTK_REAL time = 0.0; // dummy time
        ghext->amrcore->regrid(0, time);
        saved_cctkGH = nullptr;
      }
      const int new_numlevels = ghext->amrcore->finestLevel() + 1;
      const int max_numlevels = ghext->amrcore->maxLevel() + 1;
      assert(new_numlevels >= 0 && new_numlevels <= max_numlevels);
#pragma omp critical
      {
        CCTK_VINFO("  old levels %d, new levels %d", old_numlevels,
                   new_numlevels);
        double pts0 = ghext->leveldata.at(0).fab->boxArray().d_numPts();
        assert(!active_levels);
        for (const auto &leveldata : ghext->leveldata) {
          const int sz = leveldata.fab->size();
          const double pts = leveldata.fab->boxArray().d_numPts();
          if (leveldata.level == 0) {
            CCTK_VINFO("  level %d: %d boxes, %.0f cells (%.4g%%)",
                       leveldata.level, sz, pts,
                       100 * pts / (pow(2.0, dim * leveldata.level) * pts0));
          } else {
            const double ptsc = ghext->leveldata.at(leveldata.level - 1)
                                    .fab->boxArray()
                                    .d_numPts();
            CCTK_VINFO("  level %d: %d boxes, %.0f cells (%.4g%%, %.0f%%)",
                       leveldata.level, sz, pts,
                       100 * pts / (pow(2.0, dim * leveldata.level) * pts0),
                       100 * pts / (pow(2.0, dim) * ptsc));
          }
        }
      }
    }

    // Find smallest iteration number. Levels at this iteration will
    // be evolved.
    rat64 iteration = ghext->leveldata.at(0).iteration;
    for (const auto &leveldata : ghext->leveldata)
      iteration = min(iteration, leveldata.iteration);

    // Loop over all levels, in batches that combine levels that don't
    // subcycle. The level range is [min_level, max_level).
    int min_level = 0;
    while (min_level < int(ghext->leveldata.size())) {
      // Find end of batch
      int max_level = min_level + 1;
      while (max_level < int(ghext->leveldata.size()) &&
             !ghext->leveldata.at(max_level).is_subcycling_level)
        ++max_level;

      // Skip this batch of levels if it is not active at the current
      // iteration
      if (ghext->leveldata.at(min_level).iteration > iteration)
        break;

      active_levels = make_optional<active_levels_t>(min_level, max_level);

      // Advance iteration number on this batch of levels
      active_levels->loop([&](auto &restrict leveldata) {
        leveldata.iteration += leveldata.delta_iteration;
      });

      InvalidateTimelevels(cctkGH);

      CycleTimelevels(cctkGH);

      CCTK_Traverse(cctkGH, "CCTK_PRESTEP");
      CCTK_Traverse(cctkGH, "CCTK_EVOL");

      // Reflux
      assert(active_levels);
      for (int level = int(ghext->leveldata.size()) - 2; level >= 0; --level)
        Reflux(level);

      if (!restrict_during_sync) {
        // Restrict
        for (int level = int(ghext->leveldata.size()) - 2; level >= 0; --level)
          Restrict(level);
        CCTK_Traverse(cctkGH, "CCTK_POSTRESTRICT");
      }

      CCTK_Traverse(cctkGH, "CCTK_POSTSTEP");
      CCTK_Traverse(cctkGH, "CCTK_CHECKPOINT");
      CCTK_Traverse(cctkGH, "CCTK_ANALYSIS");
      CCTK_OutputGH(cctkGH);

      active_levels = optional<active_levels_t>();
    }

  } // main loop

  return 0;
}

// Schedule shutdown
int Shutdown(tFleshConfig *config) {
  assert(config);
  cGH *restrict const cctkGH = config->GH[0];
  assert(cctkGH);

  static Timer timer("Shutdown");
  Interval interval(timer);

#pragma omp critical
  CCTK_VINFO("Shutting down...");

  assert(!active_levels);
  active_levels =
      make_optional<active_levels_t>(0, int(ghext->leveldata.size()));

  CCTK_Traverse(cctkGH, "CCTK_TERMINATE");

  active_levels = optional<active_levels_t>();
  active_levels = make_optional<active_levels_t>(0, 0);

  CCTK_Traverse(cctkGH, "CCTK_SHUTDOWN");

  active_levels = optional<active_levels_t>();
  assert(!ghext);

  return 0;
}

// Call a scheduled function
int CallFunction(void *function, cFunctionData *restrict attribute,
                 void *data) {
  DECLARE_CCTK_PARAMETERS;

  assert(function);
  assert(attribute);
  assert(data);

  cGH *restrict const cctkGH = static_cast<cGH *>(data);

  if (verbose)
#pragma omp critical
    CCTK_VINFO("CallFunction iteration %d %s: %s::%s", cctkGH->cctk_iteration,
               attribute->where, attribute->thorn, attribute->routine);

  assert(active_levels);

  // Check whether input variables have valid data
  {
    const vector<clause_t> &reads = decode_clauses(attribute, rdwr_t::read);
    for (const auto &rd : reads) {
      if (CCTK_GroupTypeI(rd.gi) == CCTK_GF) {

        active_levels->loop([&](const auto &restrict leveldata) {
          const auto &restrict groupdata = *leveldata.groupdata.at(rd.gi);
          const valid_t &need = rd.valid;
          error_if_invalid(groupdata, rd.vi, rd.tl, need, [&] {
            ostringstream buf;
            buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
                << attribute->where << ": " << attribute->thorn
                << "::" << attribute->routine << " checking input";
            return buf.str();
          });
          check_valid(groupdata, rd.vi, rd.tl, [&] {
            ostringstream buf;
            buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
                << attribute->where << ": " << attribute->thorn
                << "::" << attribute->routine << " checking input";
            return buf.str();
          });
        });

      } else { // CCTK_SCALAR

        const auto &restrict scalargroupdata =
            *ghext->globaldata.scalargroupdata.at(rd.gi);
        const valid_t &need = rd.valid;
        error_if_invalid(scalargroupdata, rd.vi, rd.tl, need, [&] {
          ostringstream buf;
          buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
              << attribute->where << ": " << attribute->thorn
              << "::" << attribute->routine << " checking input";
          return buf.str();
        });
        check_valid(scalargroupdata, rd.vi, rd.tl, [&] {
          ostringstream buf;
          buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
              << attribute->where << ": " << attribute->thorn
              << "::" << attribute->routine << " checking input";
          return buf.str();
        });
      }
    }
  }

  // Poison those output variables that are not input variables
  if (poison_undefined_values) {
    map<clause_t, valid_t> isread;
    const vector<clause_t> &reads = decode_clauses(attribute, rdwr_t::read);
    for (const auto &rd : reads) {
      clause_t cl = rd;
      cl.valid = valid_t();
      assert(isread.count(cl) == 0);
      isread[cl] = rd.valid;
    }
    const vector<clause_t> &writes = decode_clauses(attribute, rdwr_t::write);
    for (const auto &wr : writes) {
      clause_t cl = wr;
      cl.valid = valid_t();
      valid_t need;
      if (isread.count(cl) > 0)
        need = isread[cl];

      if (CCTK_GroupTypeI(wr.gi) == CCTK_GF) {

        active_levels->loop([&](auto &restrict leveldata) {
          auto &restrict groupdata = *leveldata.groupdata.at(wr.gi);
          const valid_t &provided = wr.valid;
          groupdata.valid.at(wr.tl).at(wr.vi).set_and(
              need | ~provided,
              [iteration = cctkGH->cctk_iteration, where = attribute->where,
               thorn = attribute->thorn, routine = attribute->routine] {
                ostringstream buf;
                buf << "CallFunction iteration " << iteration << " " << where
                    << ": " << thorn << "::" << routine
                    << ": Poison output variables that are not input variables";
                return buf.str();
              });
          poison_invalid(groupdata, wr.vi, wr.tl);
        });

      } else { // CCTK_SCALAR

        auto &restrict scalargroupdata =
            *ghext->globaldata.scalargroupdata.at(wr.gi);
        const valid_t &provided = wr.valid;
        scalargroupdata.valid.at(wr.tl).at(wr.vi).set_and(
            need | ~provided,
            [iteration = cctkGH->cctk_iteration, where = attribute->where,
             thorn = attribute->thorn, routine = attribute->routine] {
              ostringstream buf;
              buf << "CallFunction iteration " << iteration << " " << where
                  << ": " << thorn << "::" << routine
                  << ": Poison output variables that are not input variables";
              return buf.str();
            });
        poison_invalid(scalargroupdata, wr.vi, wr.tl);
      }
    }
  }

  // Calculate checksums over variables that are not written
  checksums_t checksums;
  if (poison_undefined_values) {
    const vector<clause_t> &writes = decode_clauses(attribute, rdwr_t::write);
    const int numgroups = CCTK_NumGroups();
    vector<vector<vector<valid_t> > > gfs(numgroups);
    for (int gi = 0; gi < numgroups; ++gi) {
      const int numvars = CCTK_NumVarsInGroupI(gi);
      gfs.at(gi).resize(numvars);
      for (int vi = 0; vi < numvars; ++vi) {
        const int numtimelevels = 1; // is expanded later if necessary
        gfs.at(gi).at(vi).resize(numtimelevels);
      }
    }
    for (const auto &wr : writes) {
      if (wr.tl >= int(gfs.at(wr.gi).at(wr.vi).size()))
        gfs.at(wr.gi).at(wr.vi).resize(wr.tl + 1);
      gfs.at(wr.gi).at(wr.vi).at(wr.tl) |= wr.valid;
    }

    checksums = calculate_checksums(gfs);
  }

  const mode_t mode = decode_mode(attribute);
  switch (mode) {
  case mode_t::local: {
    // Call function once per tile
    // TODO: provide looping function
    assert(thread_local_info.empty());
    swap(thread_local_info, saved_thread_local_info);

#if 0

#pragma omp parallel
    {
      const int thread_num = omp_get_thread_num();
      thread_local_info_t &restrict thread_info =
          *thread_local_info.at(thread_num);
      cGH *restrict const threadGH = &thread_info.cctkGH;
      update_cctkGH(threadGH, cctkGH);
      TileBox &restrict thread_tilebox = thread_info.tilebox;

      // Loop over all levels
      // TODO: parallelize this loop
      active_levels->loop([&](const auto &restrict leveldata) {
        enter_level_mode(threadGH, leveldata);
        const auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling(
            {max_tile_size_x, max_tile_size_y, max_tile_size_z});
        for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid();
             ++mfi) {
          cout << "level=" << level << " mfi.currentIndex=" << mfi.tileIndex()
               << " mfi.length=" << mfi.length() << "\n";
          MFPointer mfp(mfi);
          thread_info.mfpointer = &mfp;
          enter_local_mode(threadGH, thread_tilebox, leveldata, mfp);
          CCTK_CallFunction(function, attribute, threadGH);
          leave_local_mode(threadGH, thread_tilebox, leveldata, mfp);
          thread_info.mfpointer = nullptr;
        }
        leave_level_mode(threadGH, leveldata);
      });
    }

#else

    vector<std::function<void()> > tasks;
    active_levels->loop([&](const auto &restrict leveldata) {
      const auto mfitinfo = amrex::MFItInfo().EnableTiling(
          {max_tile_size_x, max_tile_size_y, max_tile_size_z});
      // Note: The amrex::MFIter uses global variables and OpenMP barriers
      for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {
        const MFPointer mfp(mfi);

        const auto task{[level = leveldata.level, mfp, function, attribute] {
          const int thread_num = omp_get_thread_num();
          thread_local_info_t &restrict thread_info =
              *thread_local_info.at(thread_num);
          cGH *restrict const threadGH = &thread_info.cctkGH;
          TileBox &restrict thread_tilebox = thread_info.tilebox;

          const auto &restrict leveldata = ghext->leveldata.at(level);
          thread_info.mfpointer = &mfp;

          enter_level_mode(threadGH, leveldata);
          enter_local_mode(threadGH, thread_tilebox, leveldata, mfp);
          CCTK_CallFunction(function, attribute, threadGH);
          leave_local_mode(threadGH, thread_tilebox, leveldata, mfp);
          leave_level_mode(threadGH, leveldata);
          thread_info.mfpointer = nullptr;
        }};
        tasks.emplace_back(task);
      }
    });

#pragma omp parallel
    {
      // Initialize thread-local state variables
      const int thread_num = omp_get_thread_num();
      thread_local_info_t &restrict thread_info =
          *thread_local_info.at(thread_num);
      cGH *restrict const threadGH = &thread_info.cctkGH;
      update_cctkGH(threadGH, cctkGH);

      // run all tasks
#pragma omp for schedule(dynamic)
      for (size_t i = 0; i < tasks.size(); ++i)
        tasks[i]();
    }

#endif

    swap(saved_thread_local_info, thread_local_info);
    assert(thread_local_info.empty());
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

  // Check checksums
  if (poison_undefined_values)
    check_checksums(checksums);

  // Mark output variables as having valid data
  {
    const vector<clause_t> &writes = decode_clauses(attribute, rdwr_t::write);
    for (const auto &wr : writes) {
      if (CCTK_GroupTypeI(wr.gi) == CCTK_GF) {

        active_levels->loop([&](auto &restrict leveldata) {
          auto &restrict groupdata = *leveldata.groupdata.at(wr.gi);
          const valid_t &provided = wr.valid;
          groupdata.valid.at(wr.tl).at(wr.vi).set_or(
              provided,
              [iteration = cctkGH->cctk_iteration, where = attribute->where,
               thorn = attribute->thorn, routine = attribute->routine] {
                ostringstream buf;
                buf << "CallFunction iteration " << iteration << " " << where
                    << ": " << thorn << "::" << routine
                    << ": Mark output variables as valid";
                return buf.str();
              });
          check_valid(groupdata, wr.vi, wr.tl, [&]() {
            ostringstream buf;
            buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
                << attribute->where << ": " << attribute->thorn
                << "::" << attribute->routine << " checking output";
            return buf.str();
          });
        });

      } else { // CCTK_SCALAR

        auto &restrict scalargroupdata =
            *ghext->globaldata.scalargroupdata.at(wr.gi);
        const valid_t &provided = wr.valid;
        scalargroupdata.valid.at(wr.tl).at(wr.vi).set_or(
            provided,
            [iteration = cctkGH->cctk_iteration, where = attribute->where,
             thorn = attribute->thorn, routine = attribute->routine] {
              ostringstream buf;
              buf << "CallFunction iteration " << iteration << " " << where
                  << ": " << thorn << "::" << routine
                  << ": Mark output variables as valid";
              return buf.str();
            });
        check_valid(scalargroupdata, wr.vi, wr.tl, [&]() {
          ostringstream buf;
          buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
              << attribute->where << ": " << attribute->thorn
              << "::" << attribute->routine << " checking output";
          return buf.str();
        });
      }
    }
  }

  // Mark invalid variables as having invalid data
  {
    const vector<clause_t> &invalids =
        decode_clauses(attribute, rdwr_t::invalid);
    for (const auto &inv : invalids) {
      if (CCTK_GroupTypeI(inv.gi) == CCTK_GF) {

        active_levels->loop([&](auto &restrict leveldata) {
          auto &restrict groupdata = *leveldata.groupdata.at(inv.gi);
          const valid_t &provided = inv.valid;
          groupdata.valid.at(inv.tl).at(inv.vi).set_and(
              ~provided,
              [iteration = cctkGH->cctk_iteration, where = attribute->where,
               thorn = attribute->thorn, routine = attribute->routine] {
                ostringstream buf;
                buf << "CallFunction iteration " << iteration << " " << where
                    << ": " << thorn << "::" << routine
                    << ": Mark invalid variables as invalid";
                return buf.str();
              });
          check_valid(groupdata, inv.vi, inv.tl, [&]() {
            ostringstream buf;
            buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
                << attribute->where << ": " << attribute->thorn
                << "::" << attribute->routine << " checking output";
            return buf.str();
          });
        });

      } else { // CCTK_SCALAR

        auto &restrict scalargroupdata =
            *ghext->globaldata.scalargroupdata.at(inv.gi);
        const valid_t &provided = inv.valid;
        scalargroupdata.valid.at(inv.tl).at(inv.vi).set_and(
            ~provided,
            [iteration = cctkGH->cctk_iteration, where = attribute->where,
             thorn = attribute->thorn, routine = attribute->routine] {
              ostringstream buf;
              buf << "CallFunction iteration " << iteration << " " << where
                  << ": " << thorn << "::" << routine
                  << ": Mark invalid variables as invalid";
              return buf.str();
            });
        check_valid(scalargroupdata, inv.vi, inv.tl, [&]() {
          ostringstream buf;
          buf << "CallFunction iteration " << cctkGH->cctk_iteration << " "
              << attribute->where << ": " << attribute->thorn
              << "::" << attribute->routine << " checking output";
          return buf.str();
        });
      }
    }
  }

  constexpr int didsync = 0;
  return didsync;
}

int SyncGroupsByDirI(const cGH *restrict cctkGH, int numgroups,
                     const int *groups0, const int *directions) {
  DECLARE_CCTK_PARAMETERS;

  assert(in_global_mode(cctkGH));

  static Timer timer("Sync");
  Interval interval(timer);

  assert(cctkGH);
  assert(numgroups >= 0);
  assert(groups0);

  if (verbose) {
    ostringstream buf;
    for (int n = 0; n < numgroups; ++n) {
      if (n != 0)
        buf << ", ";
      buf << unique_C_ptr<char>(CCTK_GroupName(groups0[n])).get();
    }
#pragma omp critical
    CCTK_VINFO("SyncGroups %s", buf.str().c_str());
  }

  const int gi_regrid_error = CCTK_GroupIndex("CarpetX::regrid_error");
  assert(gi_regrid_error >= 0);
  const int gi_refinement_level = CCTK_GroupIndex("CarpetX::refinement_level");
  assert(gi_refinement_level >= 0);

  vector<int> groups;
  for (int n = 0; n < numgroups; ++n) {
    const int gi = groups0[n];
    if (CCTK_GroupTypeI(gi) != CCTK_GF)
      continue;
    // Don't restrict the regridding error nor the refinement level
    if (gi == gi_regrid_error || gi == gi_refinement_level)
      continue;
    groups.push_back(gi);
  }

  if (restrict_during_sync)
    active_levels->loop_reverse([&](const auto &leveldata) {
      if (leveldata.level < int(ghext->leveldata.size()) - 1)
        Restrict(leveldata.level, groups);
    });

  active_levels->loop([&](auto &restrict leveldata) {
    for (const int gi : groups) {
      cGroup group;
      int ierr = CCTK_GroupData(gi, &group);
      assert(!ierr);

      assert(group.grouptype == CCTK_GF);

      auto &restrict groupdata = *leveldata.groupdata.at(gi);
      // We always sync all directions.
      // If there is more than one time level, then we don't sync the
      // oldest.
      // TODO: during evolution, sync only one time level
      const int ntls = groupdata.mfab.size();
      const int sync_tl = ntls > 1 ? ntls - 1 : ntls;

      auto physbc_bcs = get_boundaries(groupdata);
      CarpetXPhysBCFunct &physbc = get<0>(physbc_bcs);
      const amrex::Vector<amrex::BCRec> &bcs = get<1>(physbc_bcs);

      if (leveldata.level == 0) {
        // Coarsest level: Copy from adjacent boxes on same level

        for (int tl = 0; tl < sync_tl; ++tl) {
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            // Synchronization only uses the interior
            error_if_invalid(groupdata, vi, tl, make_valid_int(),
                             [] { return "SyncGroupsByDirI before syncing"; });
            groupdata.valid.at(tl).at(vi).set_and(~make_valid_ghosts(), [] {
              return "SyncGroupsByDirI before syncing: Mark ghost zones as "
                     "invalid";
            });
            poison_invalid(groupdata, vi, tl);
            check_valid(groupdata, vi, tl,
                        [] { return "SyncGroupsByDirI before syncing"; });
          }
          {
            static Timer timer("Sync::FillPatchSingleLevel");
            Interval interval(timer);
            FillPatchSingleLevel(
                *groupdata.mfab.at(tl), 0.0, {&*groupdata.mfab.at(tl)}, {0.0},
                0, 0, groupdata.numvars, ghext->amrcore->Geom(leveldata.level),
                physbc, 0);
          }
          for (int vi = 0; vi < groupdata.numvars; ++vi)
            groupdata.valid.at(tl).at(vi).set_ghosts(true, [] {
              return "SyncGroupsByDirI after syncing: Mark ghost zones as "
                     "valid";
            });
        }

      } else {
        // Refined level: Prolongate boundaries from next coarser level, and
        // copy from adjacent boxes on same level

        const int level = leveldata.level;
        const auto &restrict coarseleveldata = ghext->leveldata.at(level - 1);
        auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
        assert(coarsegroupdata.numvars == groupdata.numvars);

        amrex::Interpolater *const interpolator =
            get_interpolator(groupdata.indextype);

        const amrex::IntVect reffact{2, 2, 2};

        for (int tl = 0; tl < sync_tl; ++tl) {
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
#warning "TODO: interpolation does not require boundaries: also for regridding!"
            error_if_invalid(coarsegroupdata, vi, tl, make_valid_int(), [] {
              return "SyncGroupsByDirI on coarse level "
                     "before prolongation";
            });
            error_if_invalid(groupdata, vi, tl, make_valid_int(), [] {
              return "SyncGroupsByDirI on fine level before prolongation";
            });
            poison_invalid(groupdata, vi, tl);
            check_valid(coarsegroupdata, vi, tl, [] {
              return "SyncGroupsByDirI on coarse level before prolongation";
            });
            check_valid(groupdata, vi, tl, [] {
              return "SyncGroupsByDirI on fine level before prolongation";
            });
            groupdata.valid.at(tl).at(vi).set_ghosts(false, [] {
              return "SyncGroupsByDirI before prolongation: Mark ghosts as "
                     "invalid";
            });
          }
          {
            static Timer timer("Sync::FillPatchTwoLevels");
            Interval interval(timer);
            FillPatchTwoLevels(
                *groupdata.mfab.at(tl), 0.0, {&*coarsegroupdata.mfab.at(tl)},
                {0.0}, {&*groupdata.mfab.at(tl)}, {0.0}, 0, 0,
                groupdata.numvars, ghext->amrcore->Geom(level - 1),
                ghext->amrcore->Geom(level), physbc, 0, physbc, 0, reffact,
                interpolator, bcs, 0);
          }
          for (int vi = 0; vi < groupdata.numvars; ++vi) {
            groupdata.valid.at(tl).at(vi).set_ghosts(
                true, [] { return "SyncGroupsByDirI after prolongation"; });
          }
        } // if not all_invalid
      }   // for tl

      for (int tl = 0; tl < sync_tl; ++tl) {
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          poison_invalid(groupdata, vi, tl);
          check_valid(groupdata, vi, tl,
                      [&]() { return "SyncGroupsByDirI after syncing"; });
        }
      }
    }
  });

  return numgroups; // number of groups synchronized
}

void Reflux(int level) {
  DECLARE_CCTK_PARAMETERS;

  if (!do_reflux)
    return;

  static Timer timer("Reflux");
  Interval interval(timer);

  auto &leveldata = ghext->leveldata.at(level);
  const auto &fineleveldata = ghext->leveldata.at(level + 1);
  for (int gi = 0; gi < int(leveldata.groupdata.size()); ++gi) {
    const int tl = 0;
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &groupdata = *leveldata.groupdata.at(gi);
    const auto &finegroupdata = *fineleveldata.groupdata.at(gi);

    // If the group has associated fluxes
    if (finegroupdata.freg) {

      // Check coarse and fine data and fluxes are valid
      for (int vi = 0; vi < finegroupdata.numvars; ++vi) {
        error_if_invalid(finegroupdata, vi, tl, make_valid_int(), [] {
          return "Reflux before refluxing: Fine level data";
        });
        error_if_invalid(groupdata, vi, tl, make_valid_int(), [] {
          return "Reflux before refluxing: Coarse level data";
        });
      }
      for (int d = 0; d < dim; ++d) {
        const int flux_gi = finegroupdata.fluxes[d];
        const auto &flux_finegroupdata = *fineleveldata.groupdata.at(flux_gi);
        const auto &flux_groupdata = *leveldata.groupdata.at(flux_gi);
        for (int vi = 0; vi < finegroupdata.numvars; ++vi) {
          error_if_invalid(flux_finegroupdata, vi, tl, make_valid_int(), [&] {
            ostringstream buf;
            buf << "Reflux: Fine level flux in direction " << d;
            return buf.str();
          });
          error_if_invalid(flux_groupdata, vi, tl, make_valid_int(), [&] {
            ostringstream buf;
            buf << "Reflux: Coarse level flux in direction " << d;
            return buf.str();
          });
        }
      }

      for (int d = 0; d < dim; ++d) {
        const int flux_gi = finegroupdata.fluxes[d];
        const auto &flux_finegroupdata = *fineleveldata.groupdata.at(flux_gi);
        const auto &flux_groupdata = *leveldata.groupdata.at(flux_gi);
        finegroupdata.freg->CrseInit(*flux_groupdata.mfab.at(tl), d, 0, 0,
                                     flux_groupdata.numvars, -1);
        finegroupdata.freg->FineAdd(*flux_finegroupdata.mfab.at(tl), d, 0, 0,
                                    flux_finegroupdata.numvars, 1);
      }
      const amrex::Geometry &geom = ghext->amrcore->Geom(level);
      finegroupdata.freg->Reflux(*groupdata.mfab.at(tl), 1.0, 0, 0,
                                 groupdata.numvars, geom);

      for (int vi = 0; vi < finegroupdata.numvars; ++vi) {
        check_valid(finegroupdata, vi, tl, [&]() {
          return "Reflux after refluxing: Fine level data";
        });
      }
    }
  } // for gi
}

void Restrict(int level, const vector<int> &groups) {
  DECLARE_CCTK_PARAMETERS;

#warning "TODO"
  assert(do_restrict);
  if (!do_restrict)
    return;

  static Timer timer("Restrict");
  Interval interval(timer);

  const int gi_regrid_error = CCTK_GroupIndex("CarpetX::regrid_error");
  assert(gi_regrid_error >= 0);
  const int gi_refinement_level = CCTK_GroupIndex("CarpetX::refinement_level");
  assert(gi_refinement_level >= 0);

  auto &leveldata = ghext->leveldata.at(level);
  const auto &fineleveldata = ghext->leveldata.at(level + 1);
  for (const int gi : groups) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    assert(group.grouptype == CCTK_GF);

    auto &groupdata = *leveldata.groupdata.at(gi);
    const auto &finegroupdata = *fineleveldata.groupdata.at(gi);
    const amrex::IntVect reffact{2, 2, 2};

    // Don't restrict the regridding error nor the refinement level
    if (gi == gi_regrid_error || gi == gi_refinement_level)
      continue;
    // Don't restrict groups that have restriction disabled
    if (!groupdata.do_restrict)
      continue;

    // If there is more than one time level, then we don't restrict the oldest.
    // TODO: during evolution, restrict only one time level
    int ntls = groupdata.mfab.size();
    int restrict_tl = ntls > 1 ? ntls - 1 : ntls;
    for (int tl = 0; tl < restrict_tl; ++tl) {

      for (int vi = 0; vi < groupdata.numvars; ++vi) {

        // Restriction only uses the interior
        error_if_invalid(finegroupdata, vi, tl, make_valid_int(), [] {
          return "Restrict on fine level before restricting";
        });
        poison_invalid(finegroupdata, vi, tl);
        check_valid(finegroupdata, vi, tl,
                    [] { return "Restrict on fine level before restricting"; });
        error_if_invalid(groupdata, vi, tl, make_valid_int(), [] {
          return "Restrict on coarse level before restricting";
        });
        poison_invalid(groupdata, vi, tl);
        check_valid(groupdata, vi, tl, [] {
          return "Restrict on coarse level before restricting";
        });
      }

#warning                                                                       \
    "TODO: Allow different restriction operators, and ensure this is conservative"
      // rank: 0: vertex, 1: edge, 2: face, 3: volume
      int rank = 0;
      for (int d = 0; d < dim; ++d)
        rank += groupdata.indextype[d];
      switch (rank) {
      case 0:
        average_down_nodal(*finegroupdata.mfab.at(tl), *groupdata.mfab.at(tl),
                           reffact);
        break;
      case 1:
        average_down_edges(*finegroupdata.mfab.at(tl), *groupdata.mfab.at(tl),
                           reffact);
        break;
      case 2:
        average_down_faces(*finegroupdata.mfab.at(tl), *groupdata.mfab.at(tl),
                           reffact);
        break;
      case 3:
        average_down(*finegroupdata.mfab.at(tl), *groupdata.mfab.at(tl), 0,
                     groupdata.numvars, reffact);
        break;
      default:
        assert(0);
      }

      // TODO: Also remember old why_valid for interior?
      for (int vi = 0; vi < groupdata.numvars; ++vi) {
        groupdata.valid.at(tl).at(vi).set(make_valid_int(),
                                          [] { return "Restrict"; });
        poison_invalid(groupdata, vi, tl);
        check_valid(groupdata, vi, tl, [&]() {
          return "Restrict on coarse level after restricting";
        });
      }

    } // for tl
  }   // for gi
}

void Restrict(int level) {
  const int numgroups = CCTK_NumGroups();
  vector<int> groups;
  groups.reserve(numgroups);
  const auto &leveldata = ghext->leveldata.at(level);
  for (const auto &groupdataptr : leveldata.groupdata) {
    // Restrict only grid functions
    if (groupdataptr) {
      auto &restrict groupdata = *groupdataptr;
      // Restrict only evolved grid functions
      if (groupdata.do_checkpoint)
        groups.push_back(groupdata.groupindex);
    }
  }
  Restrict(level, groups);
}

// storage handling
namespace {
int GroupStorageCrease(const cGH *cctkGH, int n_groups, const int *groups,
                       const int *requested_tls, int *status, const bool inc) {
  DECLARE_CCTK_PARAMETERS;

  assert(cctkGH);
  assert(n_groups >= 0);
  assert(groups);
  assert(requested_tls);
  for (int n = 0; n < n_groups; ++n) {
    if (groups[n] < 0 or groups[n] >= CCTK_NumGroups()) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Group index %d is illegal", groups[n]);
      return -1;
    }
    assert(groups[n] >= 0 and groups[n] < CCTK_NumGroups());
    assert(requested_tls[n] >= 0 or requested_tls[n] == -1);
  }

  // sanitize list of requested timelevels
  std::vector<int> tls(n_groups);
  for (int n = 0; n < n_groups; ++n) {
    int ntls = requested_tls[n];
    int const declared_tls = CCTK_DeclaredTimeLevelsGI(groups[n]);
    if (inc and declared_tls < 2 and ntls > declared_tls) {
      char *groupname = CCTK_GroupName(groups[n]);
      CCTK_VWarn(CCTK_WARN_ALERT, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "Attempting to activate %d timelevels for group '%s' which "
                 "only has a single timelevel declared in interface.ccl. "
                 "Please declared at least 2 timelevels in interface.ccl to "
                 "allow more timelevels to be created at runtime.",
                 ntls, groupname);
      free(groupname);
      ntls = declared_tls;
    }
    if (ntls == -1) {
      ntls = declared_tls;
    }
    tls.at(n) = ntls;
  }

  // TODO: actually do something
  int min_num_timelevels = INT_MAX;
  for (int n = 0; n < n_groups; ++n) {
    int const gid = groups[n];

    cGroup group;
    int ierr = CCTK_GroupData(gid, &group);
    assert(not ierr);

    // Record previous number of allocated time levels
    if (status) {
      // Note: This remembers only the last level
      status[n] = group.numtimelevels;
    }

    // Record (minimum of) current number of time levels
    min_num_timelevels = min(min_num_timelevels, group.numtimelevels);
  } // for n
  if (min_num_timelevels == INT_MAX) {
    min_num_timelevels = 0;
  }

  return min_num_timelevels;
}
} // namespace

int GroupStorageIncrease(const cGH *cctkGH, int n_groups, const int *groups,
                         const int *tls, int *status) {
  DECLARE_CCTK_PARAMETERS;

  return GroupStorageCrease(cctkGH, n_groups, groups, tls, status, true);
}

int GroupStorageDecrease(const cGH *cctkGH, int n_groups, const int *groups,
                         const int *tls, int *status) {
  DECLARE_CCTK_PARAMETERS;

  return GroupStorageCrease(cctkGH, n_groups, groups, tls, status, false);
}

int EnableGroupStorage(const cGH *cctkGH, const char *groupname) {
  const int group = CCTK_GroupIndex(groupname);
  assert(group >= 0 and group < CCTK_NumGroups());
  // TODO: decide whether to use CCTK_MaxActiveTimeLevelsGI
  const int tls = CCTK_DeclaredTimeLevelsGI(group);
  int status;
  GroupStorageIncrease(cctkGH, 1, &group, &tls, &status);
  // Return whether storage was allocated previously
  return status;
}

int DisableGroupStorage(const cGH *cctkGH, const char *groupname) {
  const int group = CCTK_GroupIndex(groupname);
  assert(group >= 0 and group < CCTK_NumGroups());
  const int tls = 0;
  int status;
  GroupStorageDecrease(cctkGH, 1, &group, &tls, &status);
  // Return whether storage was allocated previously
  return status;
}

} // namespace CarpetX
