#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <cstddef>
#include <cstdio>

namespace TestProlongate {
using namespace std;

extern "C" void TestProlongate_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    data[p.idx] = p.x;
  });
}

extern "C" void TestProlongate_Sync(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  return; // do nothing
}
extern "C" void TestProlongate_Regrid(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if (cctk_iteration < regrid_after)
    return;

  CCTK_VINFO("Setting grid at %d\n", cctk_iteration);
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    if (fabs(p.x) <= refined_radius && fabs(p.y) <= refined_radius &&
        fabs(p.z) <= refined_radius) {
      regrid_error[p.idx] = 1e3;
    } else {
      regrid_error[p.idx] = 0.;
    }
  });
}

} // namespace TestProlongate
