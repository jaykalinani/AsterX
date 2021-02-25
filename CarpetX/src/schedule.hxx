#ifndef SCHEDULE_HXX
#define SCHEDULE_HXX

#include "driver.hxx"
#include "loop.hxx"

#include <cctk.h>
#include <cctk_Schedule.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <optional>

namespace CarpetX {
using namespace std;
using namespace Loop;

int Initialise(tFleshConfig *config);
int Evolve(tFleshConfig *config);
int Shutdown(tFleshConfig *config);

int SyncGroupsByDirI(const cGH *restrict cctkGH, int numgroups,
                     const int *groups, const int *directions);

int CallFunction(void *function, cFunctionData *attribute, void *data);

int GroupStorageIncrease(const cGH *cctkGH, int n_groups, const int *groups,
                         const int *tls, int *status);
int GroupStorageDecrease(const cGH *cctkGH, int n_groups, const int *groups,
                         const int *tls, int *status);
int EnableGroupStorage(const cGH *cctkGH, const char *groupname);
int DisableGroupStorage(const cGH *cctkGH, const char *groupname);

////////////////////////////////////////////////////////////////////////////////

// This global variable passes the current cctkGH to CactusAmrCore.
// (When it is null, then CactusAmrCore does not call any scheduled
// functions. This is used early during startup.)
extern cGH *saved_cctkGH;

struct active_levels_t {
  int min_level, max_level;

  active_levels_t() = delete;

  active_levels_t(const int min_level, const int max_level)
      : min_level(min_level), max_level(max_level) {
    assert(min_level >= 0);
    assert(max_level <= int(ghext->leveldata.size()));
  }

  template <typename F> void loop(F f) {
    for (int level = min_level; level < max_level; ++level)
      assert(ghext->leveldata.at(level).iteration ==
             ghext->leveldata.at(min_level).iteration);
    for (int level = min_level; level < max_level; ++level)
      f(ghext->leveldata.at(level));
  }

  template <typename F> void loop_reverse(F f) {
    for (int level = min_level; level < max_level; ++level)
      assert(ghext->leveldata.at(level).iteration ==
             ghext->leveldata.at(min_level).iteration);
    for (int level = max_level - 1; level >= min_level; --level)
      f(ghext->leveldata.at(level));
  }
};

// The levels CallFunction should traverse
extern optional<active_levels_t> active_levels;

////////////////////////////////////////////////////////////////////////////////

// Like an MFIter, but does not support iteration, instead it can be copied
struct MFPointer {
  int m_index;
  amrex::Box m_fabbox;
  amrex::Box m_growntilebox;
  amrex::Box m_validbox;
  amrex::IntVect m_nGrowVect;

  MFPointer() = delete;
  MFPointer(const MFPointer &) = default;
  MFPointer(MFPointer &&) = default;
  MFPointer &operator=(const MFPointer &) = default;
  MFPointer &operator=(MFPointer &&) = default;
  MFPointer(const amrex::MFIter &mfi)
      : m_index((assert(mfi.isValid()), mfi.index())), m_fabbox(mfi.fabbox()),
        m_growntilebox(mfi.growntilebox()), m_validbox(mfi.validbox()),
        m_nGrowVect(mfi.theFabArrayBase().nGrowVect()) {}

  constexpr int index() const noexcept { return m_index; }
  constexpr amrex::Box fabbox() const noexcept { return m_fabbox; }
  constexpr amrex::Box growntilebox() const noexcept { return m_growntilebox; }
  constexpr amrex::Box validbox() const noexcept { return m_validbox; }
  constexpr amrex::IntVect nGrowVect() const noexcept { return m_nGrowVect; }
};

struct GridDesc : GridDescBase {

  GridDesc() = delete;
  GridDesc(const GHExt::LevelData &leveldata, const MFPointer &mfp);
  GridDesc(const cGH *cctkGH) : GridDescBase(cctkGH) {}
};

// TODO: remove this
struct GridPtrDesc : GridDesc {
  amrex::Dim3 cactus_offset;

  GridPtrDesc() = delete;
  GridPtrDesc(const GHExt::LevelData &leveldata, const MFPointer &mfp);

  template <typename T> T *ptr(const amrex::Array4<T> &vars, int vi) const {
    return vars.ptr(cactus_offset.x, cactus_offset.y, cactus_offset.z, vi);
  }
  template <typename T>
  T &idx(const amrex::Array4<T> &vars, int i, int j, int k, int vi) const {
    return vars(cactus_offset.x + i, cactus_offset.y + i, cactus_offset.z + j,
                vi);
  }
};

struct GridPtrDesc1 : GridDesc {
  amrex::Dim3 cactus_offset;
  array<int, dim> gimin, gimax;
  array<int, dim> gash;

  GridPtrDesc1() = delete;
  GridPtrDesc1(const GridPtrDesc1 &) = delete;
  GridPtrDesc1 &operator=(const GridPtrDesc1 &) = delete;

  GridPtrDesc1(const GHExt::LevelData::GroupData &groupdata,
               const MFPointer &mfp);

  template <typename T> T *ptr(const amrex::Array4<T> &vars, int vi) const {
    return vars.ptr(cactus_offset.x + gimin[0], cactus_offset.y + gimin[1],
                    cactus_offset.z + gimin[2], vi);
  }
  template <typename T>
  T &idx(const amrex::Array4<T> &vars, int i, int j, int k, int vi) const {
    return vars(cactus_offset.x + gimin[0] + i, cactus_offset.y + gimin[1] + j,
                cactus_offset.z + gimin[2] + k, vi);
  }

  template <typename T>
  GF3D1<T> gf3d(const amrex::Array4<T> &vars, int vi) const {
    return GF3D1<T>(ptr(vars, vi), gimin, gimax, gash);
  }

  friend ostream &operator<<(ostream &os, const GridPtrDesc1 &p) {
    os << "GridPtrDesc1{" << (const GridDescBase &)p << ", "
       << "cactus_offset:"
       << "{" << p.cactus_offset.x << "," << p.cactus_offset.y << ","
       << p.cactus_offset.z << "}, "
       << "gimin:"
       << "{" << p.gimin[0] << "," << p.gimin[1] << "," << p.gimin[2] << "}, "
       << "gimax:"
       << "{" << p.gimax[0] << "," << p.gimax[1] << "," << p.gimax[2] << "}, "
       << "gash:"
       << "{" << p.gash[0] << "," << p.gash[1] << "," << p.gash[2] << "}";
    return os;
  }
};

bool in_local_mode(const cGH *restrict cctkGH);
bool in_level_mode(const cGH *restrict cctkGH);
bool in_global_mode(const cGH *restrict cctkGH);
bool in_meta_mode(const cGH *restrict cctkGH);

void error_if_invalid(const GHExt::LevelData ::GroupData &grouppdata, int vi,
                      int tl, const valid_t &required,
                      const function<string()> &msg);
void warn_if_invalid(const GHExt::LevelData ::GroupData &grouppdata, int vi,
                     int tl, const valid_t &required,
                     const function<string()> &msg);
void poison_invalid(const GHExt::LevelData::GroupData &groupdata, int vi,
                    int tl);
void check_valid(const GHExt::LevelData::GroupData &groupdata, int vi, int tl,
                 const function<string()> &msg);

void error_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata,
                      int vi, int tl, const valid_t &required,
                      const function<string()> &msg);
void warn_if_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                     int tl, const valid_t &required,
                     const function<string()> &msg);
void poison_invalid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                    int tl);
void check_valid(const GHExt::GlobalData::ArrayGroupData &groupdata, int vi,
                 int tl, const function<string()> &msg);

} // namespace CarpetX

#endif // #ifndef SCHEDULE_HXX
