#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "utils.hxx"

namespace AsterX {

using namespace Arith;

extern "C" void AsterX_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  const smat<CCTK_REAL, 3, DN, DN> g {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  const CCTK_REAL detg = -1.0;
  const smat<CCTK_REAL, 3, UP, UP> invg {1.0, -3.0, 2.0, 3.0, -1.0, 0.0};

  {
    assert(calc_det(g) == detg);
    CCTK_VINFO("Test calc_det succeeded");
  }

  {
    const smat<CCTK_REAL, 3, UP, UP> invg_test = calc_inv(g, detg);
    assert(invg_test == invg);
    CCTK_VINFO("Test calc_inv succeeded");
  }
}

} // namespace AsterX
