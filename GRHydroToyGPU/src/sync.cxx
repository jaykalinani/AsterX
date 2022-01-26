#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

extern "C" void GRHydroToyGPU_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Sync;

  // do nothing
}

} // namespace GRHydroToyGPU
