#include <loop.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdio>
#include <cstdbool>
#include <cmath> 

extern "C" void SetEntropy(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTSX_SetEntropy;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Set the entropy!");

  grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          entropy(p.I) = press(p.I)/rho(p.I);

  });
}


