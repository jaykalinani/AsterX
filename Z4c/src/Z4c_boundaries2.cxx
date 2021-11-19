#include <loop_device.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_Boundaries2(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Boundaries2;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_Gamtx1(layout1, Gamtx);
  const GF3D2<CCTK_REAL> gf_Gamty1(layout1, Gamty);
  const GF3D2<CCTK_REAL> gf_Gamtz1(layout1, Gamtz);

  const Loop::GridDescBaseDevice grid(cctkGH);

  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE(const PointDesc &p)
                                    ARITH_INLINE { gf_Gamtx1(p.I) = 0; });
  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE(const PointDesc &p)
                                    ARITH_INLINE { gf_Gamty1(p.I) = 0; });
  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE(const PointDesc &p)
                                    ARITH_INLINE { gf_Gamtz1(p.I) = 0; });
}

} // namespace Z4c
