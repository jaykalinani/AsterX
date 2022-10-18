#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "utils.hxx"

namespace AsterX {

extern "C" void AsterX_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  const array<CCTK_REAL, 6> g = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  const CCTK_REAL detg = -1.0;
  const array<CCTK_REAL, 6> invg = {1.0, -3.0, 2.0, 3.0, -1.0, 0.0};

  {
    assert(calc_detg(g[0], g[1], g[2], g[3], g[4], g[5]) == detg);
    CCTK_VINFO("Test calc_detg succeeded");
  }

  {
    const array<CCTK_REAL, 6> invg_test =
        calc_upperg(g[0], g[1], g[2], g[3], g[4], g[5], detg);
    for (int i = 0; i < 6; i++)
      assert(invg_test[i] == invg[i]);
    CCTK_VINFO("Test calc_upperg succeeded");
  }
}

} // namespace AsterX
