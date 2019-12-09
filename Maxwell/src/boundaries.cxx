#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace Maxwell {

extern "C" void Maxwell_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  // do nothing
}

} // namespace Maxwell
