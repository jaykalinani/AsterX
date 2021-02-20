#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Boundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_chi1(layout1, chi);

  const GF3D2<CCTK_REAL> gf_gammatxx1(layout1, gammatxx);
  const GF3D2<CCTK_REAL> gf_gammatxy1(layout1, gammatxy);
  const GF3D2<CCTK_REAL> gf_gammatxz1(layout1, gammatxz);
  const GF3D2<CCTK_REAL> gf_gammatyy1(layout1, gammatyy);
  const GF3D2<CCTK_REAL> gf_gammatyz1(layout1, gammatyz);
  const GF3D2<CCTK_REAL> gf_gammatzz1(layout1, gammatzz);

  const GF3D2<CCTK_REAL> gf_Kh1(layout1, Kh);

  const GF3D2<CCTK_REAL> gf_Atxx1(layout1, Atxx);
  const GF3D2<CCTK_REAL> gf_Atxy1(layout1, Atxy);
  const GF3D2<CCTK_REAL> gf_Atxz1(layout1, Atxz);
  const GF3D2<CCTK_REAL> gf_Atyy1(layout1, Atyy);
  const GF3D2<CCTK_REAL> gf_Atyz1(layout1, Atyz);
  const GF3D2<CCTK_REAL> gf_Atzz1(layout1, Atzz);

  const GF3D2<CCTK_REAL> gf_Gamtx1(layout1, Gamtx);
  const GF3D2<CCTK_REAL> gf_Gamty1(layout1, Gamty);
  const GF3D2<CCTK_REAL> gf_Gamtz1(layout1, Gamtz);

  const GF3D2<CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const GF3D2<CCTK_REAL> gf_betaGx1(layout1, betaGx);
  const GF3D2<CCTK_REAL> gf_betaGy1(layout1, betaGy);
  const GF3D2<CCTK_REAL> gf_betaGz1(layout1, betaGz);

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_chi1(p.I) = 1; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatxx1(p.I) = 1; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatxy1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatxz1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatyy1(p.I) = 1; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatyz1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_gammatzz1(p.I) = 1; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Kh1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxx1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxy1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxz1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atyy1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atyz1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atzz1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Gamtx1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Gamty1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Gamtz1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Theta1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_alphaG1(p.I) = 1; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_betaGx1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_betaGy1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_betaGz1(p.I) = 0; });
}

} // namespace Z4c
