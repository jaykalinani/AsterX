#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>

#include "utils.hxx"

namespace AsterX {

using namespace Arith;

extern "C" void AsterX_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  const smat<CCTK_REAL, 3> g {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  const CCTK_REAL detg = -1.0;
  const smat<CCTK_REAL, 3> invg {1.0, -3.0, 2.0, 3.0, -1.0, 0.0};

  {
    assert(calc_det(g) == detg);
    CCTK_VINFO("Test calc_det succeeded");
  }

  {
    const smat<CCTK_REAL, 3> invg_test = calc_inv(g, detg);
    assert(invg_test(0,0) == 1.0);
    assert(invg_test(0,1) == -3.0);
    assert(invg_test(0,2) == 2.0);
    assert(invg_test(1,1) == 3.0);
    assert(invg_test(1,2) == -1.0);
    assert(invg_test(2,2) == 0.0);
    assert(invg_test.elts[0] == 1.0);
    assert(invg_test.elts[1] == -3.0);
    assert(invg_test.elts[2] == 2.0);
    assert(invg_test.elts[3] == 3.0);
    assert(invg_test.elts[4] == -1.0);
    assert(invg_test.elts[5] == 0.0);
    assert(invg_test == invg);
    CCTK_VINFO("Test calc_inv succeeded");
  }
}

} // namespace AsterX
