#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Boundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_chi_(&layout, chi);

  const GF3D2<CCTK_REAL> gf_gammatxx_(&layout, gammatxx);
  const GF3D2<CCTK_REAL> gf_gammatxy_(&layout, gammatxy);
  const GF3D2<CCTK_REAL> gf_gammatxz_(&layout, gammatxz);
  const GF3D2<CCTK_REAL> gf_gammatyy_(&layout, gammatyy);
  const GF3D2<CCTK_REAL> gf_gammatyz_(&layout, gammatyz);
  const GF3D2<CCTK_REAL> gf_gammatzz_(&layout, gammatzz);

  const GF3D2<CCTK_REAL> gf_Kh_(&layout, Kh);

  const GF3D2<CCTK_REAL> gf_Atxx_(&layout, Atxx);
  const GF3D2<CCTK_REAL> gf_Atxy_(&layout, Atxy);
  const GF3D2<CCTK_REAL> gf_Atxz_(&layout, Atxz);
  const GF3D2<CCTK_REAL> gf_Atyy_(&layout, Atyy);
  const GF3D2<CCTK_REAL> gf_Atyz_(&layout, Atyz);
  const GF3D2<CCTK_REAL> gf_Atzz_(&layout, Atzz);

  const GF3D2<CCTK_REAL> gf_Gamtx_(&layout, Gamtx);
  const GF3D2<CCTK_REAL> gf_Gamty_(&layout, Gamty);
  const GF3D2<CCTK_REAL> gf_Gamtz_(&layout, Gamtz);

  const GF3D2<CCTK_REAL> gf_Theta_(&layout, Theta);

  const GF3D2<CCTK_REAL> gf_alphaG_(&layout, alphaG);

  const GF3D2<CCTK_REAL> gf_betaGx_(&layout, betaGx);
  const GF3D2<CCTK_REAL> gf_betaGy_(&layout, betaGy);
  const GF3D2<CCTK_REAL> gf_betaGz_(&layout, betaGz);

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_chi_(p.I) = 1; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatxx_(p.I) = 1; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatxy_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatxz_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatyy_(p.I) = 1; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatyz_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatzz_(p.I) = 1; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Kh_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxx_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxy_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxz_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atyy_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atyz_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atzz_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Gamtx_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Gamty_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Gamtz_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Theta_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_alphaG_(p.I) = 1; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_betaGx_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_betaGy_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_betaGz_(p.I) = 0; });
}

} // namespace Z4c
