#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "utils.hxx"

namespace AsterX {

extern "C" void Test_utils(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  {
    array<CCTK_REAL, 6> g = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    assert(calc_detg(g[0], g[1], g[2], g[3], g[4], g[5]) == -1.0);
    CCTK_VINFO("Test calc_detg succeeded");
  }
}

} // namespace AsterX
