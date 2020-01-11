#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Boundaries;

  const GF3D<CCTK_REAL, 0, 0, 0> gf_chi_(cctkGH, chi);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxx_(cctkGH, gammatxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxy_(cctkGH, gammatxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxz_(cctkGH, gammatxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyy_(cctkGH, gammatyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyz_(cctkGH, gammatyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatzz_(cctkGH, gammatzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Kh_(cctkGH, Kh);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxx_(cctkGH, Atxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxy_(cctkGH, Atxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxz_(cctkGH, Atxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyy_(cctkGH, Atyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyz_(cctkGH, Atyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atzz_(cctkGH, Atzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamtx_(cctkGH, Gamtx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamty_(cctkGH, Gamty);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamtz_(cctkGH, Gamtz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Theta_(cctkGH, Theta);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_alphaG_(cctkGH, alphaG);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGx_(cctkGH, betaGx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGy_(cctkGH, betaGy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGz_(cctkGH, betaGz);

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
