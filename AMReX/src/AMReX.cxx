#include <AMReX.hxx>

#include <AMReX.H>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <mpi.h>

#include <string>
#include <type_traits>
#include <utility>

namespace AMReX {
using namespace std;

amrex::AMReX *pamrex = nullptr;
unique_ptr<GHExt> ghext;

extern "C" int AMReX_Startup() {
  string banner = "AMR driver provided by AMReX " + amrex::Version();
  CCTK_RegisterBanner(banner.c_str());

  // Initialize AMReX
  pamrex = amrex::Initialize(MPI_COMM_WORLD);

  // Create grid structure
  ghext = make_unique<GHExt>();

  // Define box array
  IntVect dom_lo(AMREX_D_DECL(0, 0, 0));
  IntVect dom_hi(
      AMREX_D_DECL(ghext->ncells - 1, ghext->ncells - 1, ghext->ncells - 1));
  Box domain(dom_lo, dom_hi);
  ghext->ba.define(domain);

  // Break up box array into chunks no larger than max_grid_size along
  // each direction
  const int max_grid_size = 32;
  ghext->ba.maxSize(max_grid_size);

  // Define physical box
  RealBox real_box({AMREX_D_DECL(-1.0, -1.0, -1.0)},
                   {AMREX_D_DECL(1.0, 1.0, 1.0)});

  // Define geometry
  Vector<int> is_periodic(AMREX_SPACEDIM, 1); // periodic in all directions
  ghext->geom.define(domain, &real_box, CoordSys::cartesian,
                     is_periodic.data());

  const int nvars = 4; // number of grid functions

  // Distributed boxes
  DistributionMapping dm(ghext->ba);

  // Allocate grid hierarchy
  ghext->mfab = MultiFab(ghext->ba, dm, nvars, ghext->nghostzones);

  return 0;
}

extern "C" int AMReX_Shutdown() {
  // Deallocate grid hierarchy
  ghext = nullptr;

  // Finalize AMReX
  amrex::Finalize(pamrex);
  pamrex = nullptr;

  return 0;
}

} // namespace AMReX
