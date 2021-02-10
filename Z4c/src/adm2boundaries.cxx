#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_ADM2Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM2Boundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};

  const GF3D1<CCTK_REAL> gf_dtkxx_(cctkGH, indextype, nghostzones, dtkxx);
  const GF3D1<CCTK_REAL> gf_dtkxy_(cctkGH, indextype, nghostzones, dtkxy);
  const GF3D1<CCTK_REAL> gf_dtkxz_(cctkGH, indextype, nghostzones, dtkxz);
  const GF3D1<CCTK_REAL> gf_dtkyy_(cctkGH, indextype, nghostzones, dtkyy);
  const GF3D1<CCTK_REAL> gf_dtkyz_(cctkGH, indextype, nghostzones, dtkyz);
  const GF3D1<CCTK_REAL> gf_dtkzz_(cctkGH, indextype, nghostzones, dtkzz);

  const GF3D1<CCTK_REAL> gf_dt2alp_(cctkGH, indextype, nghostzones, dt2alp);

  const GF3D1<CCTK_REAL> gf_dt2betax_(cctkGH, indextype, nghostzones, dt2betax);
  const GF3D1<CCTK_REAL> gf_dt2betay_(cctkGH, indextype, nghostzones, dt2betay);
  const GF3D1<CCTK_REAL> gf_dt2betaz_(cctkGH, indextype, nghostzones, dt2betaz);

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
