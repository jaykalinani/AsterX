#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cassert>

namespace HydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

extern "C" void HydroToyGPU_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_Boundaries;

  // Boundary conditions are not implemented yet. We require a grid
  // structure that has no boundaries, i.e. which has symmetries on
  // all boundaries.

  // do nothing

  const GridDescBaseDevice grid(cctkGH);

  grid.loop_bnd_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE {
                              assert(false); // This should not be executed
                            });
}

} // namespace HydroToyGPU
