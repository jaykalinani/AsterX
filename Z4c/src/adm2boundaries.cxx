#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_ADM2Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM2Boundaries;

  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkxx_(cctkGH, dtkxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkxy_(cctkGH, dtkxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkxz_(cctkGH, dtkxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkyy_(cctkGH, dtkyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkyz_(cctkGH, dtkyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtkzz_(cctkGH, dtkzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_dt2alp_(cctkGH, dt2alp);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_dt2betax_(cctkGH, dt2betax);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dt2betay_(cctkGH, dt2betay);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dt2betaz_(cctkGH, dt2betaz);

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
