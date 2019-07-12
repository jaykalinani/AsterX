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

namespace AMReX {
using namespace std;

static_assert(is_same<Real, CCTK_REAL>::value,
              "AMReX's Real type must be the same as Cactus's CCTK_REAL");

struct GHExt {
  const int ncells = 100;    // grid size
  const int nghostzones = 1; // number of ghost zones

  CCTK_REAL time, delta_time;

  BoxArray ba;
  Geometry geom;
  MultiFab mfab;
};

extern unique_ptr<GHExt> ghext;

} // namespace AMReX

#endif // #ifndef AMREX_HXX
