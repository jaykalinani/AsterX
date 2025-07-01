#include <loop.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

extern "C" void SetYe(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTSX_SetYe;
  DECLARE_CCTK_PARAMETERS;

  if (set_Ye_postinitial) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          Ye(p.I) = Ye_atmosphere;
        });
  }

}
