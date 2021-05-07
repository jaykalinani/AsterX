#include <loop_device.hxx>

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_ConstraintBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ConstraintBoundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_ZtCx1(layout1, ZtCx);
  const GF3D2<CCTK_REAL> gf_ZtCy1(layout1, ZtCy);
  const GF3D2<CCTK_REAL> gf_ZtCz1(layout1, ZtCz);

  const GF3D2<CCTK_REAL> gf_HC1(layout1, HC);

  const GF3D2<CCTK_REAL> gf_MtCx1(layout1, MtCx);
  const GF3D2<CCTK_REAL> gf_MtCy1(layout1, MtCy);
  const GF3D2<CCTK_REAL> gf_MtCz1(layout1, MtCz);

  const GF3D2<CCTK_REAL> gf_allC1(layout1, allC);

  const Loop::GridDescBaseDevice grid(cctkGH);

  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p)
                                    ARITH_INLINE { gf_ZtCx1(p.I) = 0; });
  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p)
                                    ARITH_INLINE { gf_ZtCy1(p.I) = 0; });
  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p)
                                    ARITH_INLINE { gf_ZtCz1(p.I) = 0; });

  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p)
                                    ARITH_INLINE { gf_HC1(p.I) = 0; });

  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p)
                                    ARITH_INLINE { gf_MtCx1(p.I) = 0; });
  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p)
                                    ARITH_INLINE { gf_MtCy1(p.I) = 0; });
  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p)
                                    ARITH_INLINE { gf_MtCz1(p.I) = 0; });

  grid.loop_bnd_device<0, 0, 0>(grid.nghostzones,
                                [=] ARITH_DEVICE ARITH_HOST(const PointDesc &p)
                                    ARITH_INLINE { gf_allC1(p.I) = 0; });
}

} // namespace Z4c
