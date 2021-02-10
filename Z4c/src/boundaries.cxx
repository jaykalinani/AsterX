#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Boundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};

  const GF3D1<CCTK_REAL> gf_chi_(cctkGH, indextype, nghostzones, chi);

  const GF3D1<CCTK_REAL> gf_gammatxx_(cctkGH, indextype, nghostzones, gammatxx);
  const GF3D1<CCTK_REAL> gf_gammatxy_(cctkGH, indextype, nghostzones, gammatxy);
  const GF3D1<CCTK_REAL> gf_gammatxz_(cctkGH, indextype, nghostzones, gammatxz);
  const GF3D1<CCTK_REAL> gf_gammatyy_(cctkGH, indextype, nghostzones, gammatyy);
  const GF3D1<CCTK_REAL> gf_gammatyz_(cctkGH, indextype, nghostzones, gammatyz);
  const GF3D1<CCTK_REAL> gf_gammatzz_(cctkGH, indextype, nghostzones, gammatzz);

  const GF3D1<CCTK_REAL> gf_Kh_(cctkGH, indextype, nghostzones, Kh);

  const GF3D1<CCTK_REAL> gf_Atxx_(cctkGH, indextype, nghostzones, Atxx);
  const GF3D1<CCTK_REAL> gf_Atxy_(cctkGH, indextype, nghostzones, Atxy);
  const GF3D1<CCTK_REAL> gf_Atxz_(cctkGH, indextype, nghostzones, Atxz);
  const GF3D1<CCTK_REAL> gf_Atyy_(cctkGH, indextype, nghostzones, Atyy);
  const GF3D1<CCTK_REAL> gf_Atyz_(cctkGH, indextype, nghostzones, Atyz);
  const GF3D1<CCTK_REAL> gf_Atzz_(cctkGH, indextype, nghostzones, Atzz);

  const GF3D1<CCTK_REAL> gf_Gamtx_(cctkGH, indextype, nghostzones, Gamtx);
  const GF3D1<CCTK_REAL> gf_Gamty_(cctkGH, indextype, nghostzones, Gamty);
  const GF3D1<CCTK_REAL> gf_Gamtz_(cctkGH, indextype, nghostzones, Gamtz);

  const GF3D1<CCTK_REAL> gf_Theta_(cctkGH, indextype, nghostzones, Theta);

  const GF3D1<CCTK_REAL> gf_alphaG_(cctkGH, indextype, nghostzones, alphaG);

  const GF3D1<CCTK_REAL> gf_betaGx_(cctkGH, indextype, nghostzones, betaGx);
  const GF3D1<CCTK_REAL> gf_betaGy_(cctkGH, indextype, nghostzones, betaGy);
  const GF3D1<CCTK_REAL> gf_betaGz_(cctkGH, indextype, nghostzones, betaGz);

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
