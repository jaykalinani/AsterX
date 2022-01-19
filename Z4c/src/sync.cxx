#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments.h>

namespace Z4c {

extern "C" void Z4c_Sync(CCTK_ARGUMENTS) { DECLARE_CCTK_ARGUMENTS_Z4c_Sync; }

} // namespace Z4c
