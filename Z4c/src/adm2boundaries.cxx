#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_ADM2Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM2Boundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_dtkxx_(&layout, dtkxx);
  const GF3D2<CCTK_REAL> gf_dtkxy_(&layout, dtkxy);
  const GF3D2<CCTK_REAL> gf_dtkxz_(&layout, dtkxz);
  const GF3D2<CCTK_REAL> gf_dtkyy_(&layout, dtkyy);
  const GF3D2<CCTK_REAL> gf_dtkyz_(&layout, dtkyz);
  const GF3D2<CCTK_REAL> gf_dtkzz_(&layout, dtkzz);

  const GF3D2<CCTK_REAL> gf_dt2alp_(&layout, dt2alp);

  const GF3D2<CCTK_REAL> gf_dt2betax_(&layout, dt2betax);
  const GF3D2<CCTK_REAL> gf_dt2betay_(&layout, dt2betay);
  const GF3D2<CCTK_REAL> gf_dt2betaz_(&layout, dt2betaz);

  //

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkxx_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkxy_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkxz_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkyy_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkyz_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtkzz_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dt2alp_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dt2betax_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dt2betay_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dt2betaz_(p.I) = 0; });
}

} // namespace Z4c
