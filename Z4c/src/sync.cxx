#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_Sync(CCTK_ARGUMENTS) { DECLARE_CCTK_ARGUMENTS_Z4c_Sync; }

} // namespace Z4c
