#ifndef SCHEDULE_HXX
#define SCHEDULE_HXX

#include "driver.hxx"
#include "loop.hxx"

#include <cctk.h>
#include <cctk_Schedule.h>
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#include <algorithm>
#include <array>

namespace AMReX {
using namespace Loop;
using namespace std;
using std::max;
using std::min;

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

// Whether CallFunction traverses all levels (-1) or just one specific level
// (>=0)
extern int current_level;

////////////////////////////////////////////////////////////////////////////////

struct GridDesc : GridDescBase {

  GridDesc() = delete;
  GridDesc(const GHExt::LevelData &leveldata, const MFIter &mfi);
  GridDesc(const cGH *cctkGH) : GridDescBase(cctkGH) {}
};

// TODO: remove this
struct GridPtrDesc : GridDesc {
  Dim3 cactus_offset;

  GridPtrDesc() = delete;
  GridPtrDesc(const GHExt::LevelData &leveldata, const MFIter &mfi);

  template <typename T> T *ptr(const Array4<T> &vars, int vi) const {
    return vars.ptr(cactus_offset.x, cactus_offset.y, cactus_offset.z, vi);
  }
  template <typename T>
  T &idx(const Array4<T> &vars, int i, int j, int k, int vi) const {
    return vars(cactus_offset.x + i, cactus_offset.y + i, cactus_offset.z + j,
                vi);
  }
};

struct GridPtrDesc1 : GridDesc {
  Dim3 cactus_offset;
  array<int, dim> gimin, gimax;
  array<int, dim> gash;

  GridPtrDesc1() = delete;
  GridPtrDesc1(const GHExt::LevelData &leveldata,
               const GHExt::LevelData::GroupData &groupdata, const MFIter &mfi);

  template <typename T> T *ptr(const Array4<T> &vars, int vi) const {
    return vars.ptr(cactus_offset.x + gimin[0], cactus_offset.y + gimin[1],
                    cactus_offset.z + gimin[2], vi);
  }
  template <typename T>
  T &idx(const Array4<T> &vars, int i, int j, int k, int vi) const {
    return vars(cactus_offset.x + gimin[0] + i, cactus_offset.y + gimin[1] + i,
                cactus_offset.z + gimin[2] + j, vi);
  }

  template <typename T> GF3D1<T> gf3d(const Array4<T> &vars, int vi) const {
    return GF3D1<T>(ptr(vars, vi), gimin, gimax, gash);
  }
};

void poison_invalid(const GHExt::LevelData &leveldata,
                    const GHExt::LevelData::GroupData &groupdata, int vi,
                    int tl);
void check_valid(const GHExt::LevelData &leveldata,
                 const GHExt::LevelData::GroupData &groupdata, int vi, int tl,
                 const function<string()> &msg);

} // namespace AMReX

#endif // #ifndef SCHEDULE_HXX
