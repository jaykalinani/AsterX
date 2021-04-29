#include "driver.hxx"
#include "loop.hxx"
#include "schedule.hxx"
#include "valid.hxx"

#include <cctk.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <functional>
#include <limits>
#include <mutex>
#include <string>

namespace CarpetX {
using namespace std;

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
  const auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
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
  const auto nan_update = [&](const GridDescBase &grid,
                              const Loop::PointDesc &p) {
    // #pragma omp critical
    static mutex m;
    lock_guard<mutex> g(m);
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
  };
  const auto nan_check = [&](const GridDescBase &grid,
                             const GF3D1<const CCTK_REAL> &ptr_,
                             const Loop::PointDesc &p) {
    if (CCTK_BUILTIN_EXPECT(!CCTK_isfinite(ptr_(p.I)), false))
      nan_update(grid, p);
  };
  const auto &leveldata = groupdata.leveldata();
  const auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
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
          "%s: Grid function \"%s\" has %td nans or infinities on refinement "
          "level %d, time level %d, in box [%d,%d,%d]:[%d,%d,%d] "
          "(%g,%g,%g):(%g,%g,%g); expected valid %s",
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

      CCTK_VERROR("%s: Grid function \"%s\" has nans or infinities on "
                  "refinement level %d, time level %d; expected valid %s",
                  msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
                  leveldata.level, tl,
                  string(groupdata.valid.at(tl).at(vi)).c_str());
    }
  }
}

// Ensure arrays are valid
void error_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata,
                      int vi, int tl, const valid_t &required,
                      const function<string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VERROR("%s: Array \"%s\" is invalid on time level %d; "
                "required %s, found %s",
                msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
                tl, string(required).c_str(),
                string(groupdata.valid.at(tl).at(vi)).c_str());
}
void warn_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                     int tl, const valid_t &required,
                     const function<string()> &msg) {
  const valid_t &have = groupdata.valid.at(tl).at(vi).get();
  if (CCTK_BUILTIN_EXPECT((required & ~have).valid_any(), false))
    CCTK_VWARN(CCTK_WARN_ALERT,
               "%s: Array \"%s\" is invalid on time level %d; "
               "required %s, found %s",
               msg().c_str(), CCTK_FullVarName(groupdata.firstvarindex + vi),
               tl, string(required).c_str(),
               string(groupdata.valid.at(tl).at(vi)).c_str());
}

// Set arrays to nan
void poison_invalid(const GHExt::GlobalData::ArrayGroupData &arraygroupdata,
                    int vi, int tl) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = arraygroupdata.valid.at(tl).at(vi).get();
  if (valid.valid_all())
    return;

  if (!valid.valid_int) {
    int dimension = arraygroupdata.dimension;
    CCTK_REAL *restrict const ptr =
        const_cast<CCTK_REAL *>(&arraygroupdata.data.at(tl).at(vi));
    const int *gsh = arraygroupdata.gsh;
    int n_elems = 1;
    for (int i = 0; i < dimension; i++)
      n_elems *= gsh[i];
    for (int i = 0; i < n_elems; i++)
      ptr[i] = 0.0 / 0.0;
  }
}

// Ensure arrays are not nan
void check_valid(const GHExt::GlobalData::ArrayGroupData &arraygroupdata,
                 int vi, int tl, const function<string()> &msg) {
  DECLARE_CCTK_PARAMETERS;
  if (!poison_undefined_values)
    return;

  const valid_t &valid = arraygroupdata.valid.at(tl).at(vi).get();
  if (!valid.valid_any())
    return;

  // arrays have no boundary so we expect them to alway be valid
  assert(valid.valid_outer && valid.valid_ghosts);

  atomic<size_t> nan_count{0};
  if (valid.valid_int) {
    const CCTK_REAL *restrict const ptr = &arraygroupdata.data.at(tl).at(vi);
    int dimension = arraygroupdata.dimension;
    const int *gsh = arraygroupdata.gsh;
    int n_elems = 1;
    for (int i = 0; i < dimension; i++)
      n_elems *= gsh[i];
    for (int i = 0; i < n_elems; i++) {
      if (CCTK_BUILTIN_EXPECT(!CCTK_isfinite(ptr[i]), false)) {
        ++nan_count;
      }
    }
  }

  if (CCTK_BUILTIN_EXPECT(nan_count > 0, false))
    CCTK_VERROR("%s: Array \"%s\" has %td nans on time level %d; "
                "expected valid %s",
                msg().c_str(),
                CCTK_FullVarName(arraygroupdata.firstvarindex + vi),
                size_t(nan_count), tl,
                string(arraygroupdata.valid.at(tl).at(vi)).c_str());
}

checksums_t
calculate_checksums(const vector<vector<vector<valid_t> > > &will_write) {
  DECLARE_CCTK_PARAMETERS;

  checksums_t checksums;

  if (!poison_undefined_values)
    return checksums;

  assert(active_levels);
  active_levels->loop([&](auto &restrict leveldata) {
    auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
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

void check_checksums(const checksums_t &checksums) {
  DECLARE_CCTK_PARAMETERS;

  if (!poison_undefined_values)
    return;
  if (checksums.empty())
    return;

  assert(active_levels);
  active_levels->loop([&](auto &restrict leveldata) {
    auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
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

} // namespace CarpetX
