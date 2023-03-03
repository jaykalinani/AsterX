#include "defs.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Hydro {
using namespace std;

extern "C" void Hydro_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Hydro_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  // do nothing
}

} // namespace Hydro
