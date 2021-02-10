#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHSBoundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};

  const GF3D1<CCTK_REAL> gf_chi_rhs_(cctkGH, indextype, nghostzones, chi_rhs);

  const GF3D1<CCTK_REAL> gf_gammatxx_rhs_(cctkGH, indextype, nghostzones,
                                          gammatxx_rhs);
  const GF3D1<CCTK_REAL> gf_gammatxy_rhs_(cctkGH, indextype, nghostzones,
                                          gammatxy_rhs);
  const GF3D1<CCTK_REAL> gf_gammatxz_rhs_(cctkGH, indextype, nghostzones,
                                          gammatxz_rhs);
  const GF3D1<CCTK_REAL> gf_gammatyy_rhs_(cctkGH, indextype, nghostzones,
                                          gammatyy_rhs);
  const GF3D1<CCTK_REAL> gf_gammatyz_rhs_(cctkGH, indextype, nghostzones,
                                          gammatyz_rhs);
  const GF3D1<CCTK_REAL> gf_gammatzz_rhs_(cctkGH, indextype, nghostzones,
                                          gammatzz_rhs);

  const GF3D1<CCTK_REAL> gf_Kh_rhs_(cctkGH, indextype, nghostzones, Kh_rhs);

  const GF3D1<CCTK_REAL> gf_Atxx_rhs_(cctkGH, indextype, nghostzones, Atxx_rhs);
  const GF3D1<CCTK_REAL> gf_Atxy_rhs_(cctkGH, indextype, nghostzones, Atxy_rhs);
  const GF3D1<CCTK_REAL> gf_Atxz_rhs_(cctkGH, indextype, nghostzones, Atxz_rhs);
  const GF3D1<CCTK_REAL> gf_Atyy_rhs_(cctkGH, indextype, nghostzones, Atyy_rhs);
  const GF3D1<CCTK_REAL> gf_Atyz_rhs_(cctkGH, indextype, nghostzones, Atyz_rhs);
  const GF3D1<CCTK_REAL> gf_Atzz_rhs_(cctkGH, indextype, nghostzones, Atzz_rhs);

  const GF3D1<CCTK_REAL> gf_Gamtx_rhs_(cctkGH, indextype, nghostzones,
                                       Gamtx_rhs);
  const GF3D1<CCTK_REAL> gf_Gamty_rhs_(cctkGH, indextype, nghostzones,
                                       Gamty_rhs);
  const GF3D1<CCTK_REAL> gf_Gamtz_rhs_(cctkGH, indextype, nghostzones,
                                       Gamtz_rhs);

  const GF3D1<CCTK_REAL> gf_Theta_rhs_(cctkGH, indextype, nghostzones,
                                       Theta_rhs);

  const GF3D1<CCTK_REAL> gf_alphaG_rhs_(cctkGH, indextype, nghostzones,
                                        alphaG_rhs);

  const GF3D1<CCTK_REAL> gf_betaGx_rhs_(cctkGH, indextype, nghostzones,
                                        betaGx_rhs);
  const GF3D1<CCTK_REAL> gf_betaGy_rhs_(cctkGH, indextype, nghostzones,
                                        betaGy_rhs);
  const GF3D1<CCTK_REAL> gf_betaGz_rhs_(cctkGH, indextype, nghostzones,
                                        betaGz_rhs);

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_chi_rhs_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatxx_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatxy_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatxz_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatyy_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatyz_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatzz_rhs_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Kh_rhs_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxx_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxy_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxz_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atyy_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atyz_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atzz_rhs_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_Gamtx_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_Gamty_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_Gamtz_rhs_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_Theta_rhs_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_alphaG_rhs_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_betaGx_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_betaGy_rhs_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_betaGz_rhs_(p.I) = 0; });
}

} // namespace Z4c
