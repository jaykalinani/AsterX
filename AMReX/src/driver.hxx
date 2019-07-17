#ifndef DRIVER_HXX
#define DRIVER_HXX

#include <AMReX.H>
#include <AMReX_AmrMesh.H>
#include <AMReX_MultiFab.H>

#include <cctk.h>

#include <memory>
#include <type_traits>
#include <vector>

namespace AMReX {
using namespace amrex;
using namespace std;

constexpr int dim = 3;

static_assert(AMREX_SPACEDIM == dim,
              "AMReX's AMREX_SPACEDIM must be the same as Cactus's cctk_dim");

static_assert(is_same<Real, CCTK_REAL>::value,
              "AMReX's Real type must be the same as Cactus's CCTK_REAL");

struct GHExt {

  // AMReX grid structure
  unique_ptr<AmrMesh> amrmesh;

  struct LevelData {
    int level;

    struct GroupData {
      int firstvarindex;
      int numvars;
      // each MultiFab has numvars components
      vector<unique_ptr<MultiFab> > mfab; // [time level]
    };
    vector<GroupData> groupdata;
  };
  vector<LevelData> leveldata;
};

extern unique_ptr<GHExt> ghext;

} // namespace AMReX

#endif // #ifndef DRIVER_HXX
