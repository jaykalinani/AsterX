#ifndef SCHEDULE_HXX
#define SCHEDULE_HXX

#include "driver.hxx"
#include "loop.hxx"

#include <cctk.h>
#include <cctk_Schedule.h>

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

////////////////////////////////////////////////////////////////////////////////

struct GridDesc : GridDescBase {

  GridDesc() = delete;
  GridDesc(const GHExt::LevelData &leveldata, const MFIter &mfi);
  GridDesc(const cGH *cctkGH) : GridDescBase(cctkGH) {}
};

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

} // namespace AMReX

#endif // #ifndef SCHEDULE_HXX
