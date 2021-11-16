#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cassert>

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

extern "C" void GRHydroToyGPU_FluxBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_FluxBoundaries;

  // Flux boundary conditions are not implemented yet. We require a
  // grid structure that has no boundaries, i.e. which has symmetries
  // on all boundaries.

  // do nothing

  const GridDescBaseDevice grid(cctkGH);

  grid.loop_bnd_device<0, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                              assert(false); // This should not be executed
                            });
  grid.loop_bnd_device<1, 0, 1>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                              assert(false); // This should not be executed
                            });
  grid.loop_bnd_device<1, 1, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                              assert(false); // This should not be executed
                            });
}

} // namespace GRHydroToyGPU
