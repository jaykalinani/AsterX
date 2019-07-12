#ifndef AMREX_HXX
#define AMREX_HXX

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
using namespace amrex;

#include <cctk.h>

#include <memory>
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
  // Parameters that should not be here
  const int ncells = 100;    // grid size
  const int nghostzones = 1; // number of ghost zones

  // AMReX grid structure
  BoxArray ba;
  Geometry geom;
  MultiFab mfab;
};

extern unique_ptr<GHExt> ghext;

// This needs to go away again!
extern vector<MFIter *> mfis;

} // namespace AMReX

#endif // #ifndef AMREX_HXX
