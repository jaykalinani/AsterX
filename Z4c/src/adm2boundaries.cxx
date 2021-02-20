#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_ADM2Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM2Boundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_dtkxx1(layout1, dtkxx);
  const GF3D2<CCTK_REAL> gf_dtkxy1(layout1, dtkxy);
  const GF3D2<CCTK_REAL> gf_dtkxz1(layout1, dtkxz);
  const GF3D2<CCTK_REAL> gf_dtkyy1(layout1, dtkyy);
  const GF3D2<CCTK_REAL> gf_dtkyz1(layout1, dtkyz);
  const GF3D2<CCTK_REAL> gf_dtkzz1(layout1, dtkzz);

  const GF3D2<CCTK_REAL> gf_dt2alp1(layout1, dt2alp);

  const GF3D2<CCTK_REAL> gf_dt2betax1(layout1, dt2betax);
  const GF3D2<CCTK_REAL> gf_dt2betay1(layout1, dt2betay);
  const GF3D2<CCTK_REAL> gf_dt2betaz1(layout1, dt2betaz);

  //

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkxx1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkxy1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkxz1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkyy1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkyz1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkzz1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dt2alp1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dt2betax1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dt2betay1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dt2betaz1(p.I) = 0; });
}

} // namespace Z4c
