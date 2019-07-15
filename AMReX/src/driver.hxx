#ifndef DRIVER_HXX
#define DRIVER_HXX

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
using namespace amrex;

#include <cctk.h>

#include <memory>
#include <type_traits>
#include <vector>

namespace AMReX {
using namespace std;

constexpr int dim = 3;

static_assert(
    AMREX_SPACEDIM == dim,
    "AMReX's number of dimensions must be the same as Cactus's cctk_dim");

static_assert(is_same<Real, CCTK_REAL>::value,
              "AMReX's Real type must be the same as Cactus's CCTK_REAL");

struct GHExt {

  struct LevelData {
    int level;

    // AMReX grid structure
    BoxArray grids;
    Geometry geom;
    DistributionMapping dmap;

    struct GroupData {
      int firstvarindex;
      int numvars;
      vector<unique_ptr<MultiFab> > mfab; // [time level]
    };
    vector<GroupData> groupdata;
  };
  vector<LevelData> leveldata;
};

extern unique_ptr<GHExt> ghext;

} // namespace AMReX

#endif // #ifndef DRIVER_HXX
