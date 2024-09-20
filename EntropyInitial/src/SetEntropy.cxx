#include <loop.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cstdio>
#include <cstdbool>
#include <cmath> 

#include "eos.hxx"
#include "eos_idealgas.hxx"

extern "C" void SetEntropy(CCTK_ARGUMENTS)
{

  using namespace EOSX;

  DECLARE_CCTK_ARGUMENTSX_SetEntropy;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Set the entropy!");

  grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

          entropy(p.I) = eos_th.entropy_from_valid_rho_eps_ye(rho(p.I),eps(p.I),1.0);

  });
}


