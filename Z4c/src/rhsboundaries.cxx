#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHSBoundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_chi_rhs_(&layout, chi_rhs);

  const GF3D2<CCTK_REAL> gf_gammatxx_rhs_(&layout,
                                          gammatxx_rhs);
  const GF3D2<CCTK_REAL> gf_gammatxy_rhs_(&layout,
                                          gammatxy_rhs);
  const GF3D2<CCTK_REAL> gf_gammatxz_rhs_(&layout,
                                          gammatxz_rhs);
  const GF3D2<CCTK_REAL> gf_gammatyy_rhs_(&layout,
                                          gammatyy_rhs);
  const GF3D2<CCTK_REAL> gf_gammatyz_rhs_(&layout,
                                          gammatyz_rhs);
  const GF3D2<CCTK_REAL> gf_gammatzz_rhs_(&layout,
                                          gammatzz_rhs);

  const GF3D2<CCTK_REAL> gf_Kh_rhs_(&layout, Kh_rhs);

  const GF3D2<CCTK_REAL> gf_Atxx_rhs_(&layout, Atxx_rhs);
  const GF3D2<CCTK_REAL> gf_Atxy_rhs_(&layout, Atxy_rhs);
  const GF3D2<CCTK_REAL> gf_Atxz_rhs_(&layout, Atxz_rhs);
  const GF3D2<CCTK_REAL> gf_Atyy_rhs_(&layout, Atyy_rhs);
  const GF3D2<CCTK_REAL> gf_Atyz_rhs_(&layout, Atyz_rhs);
  const GF3D2<CCTK_REAL> gf_Atzz_rhs_(&layout, Atzz_rhs);

  const GF3D2<CCTK_REAL> gf_Gamtx_rhs_(&layout,
                                       Gamtx_rhs);
  const GF3D2<CCTK_REAL> gf_Gamty_rhs_(&layout,
                                       Gamty_rhs);
  const GF3D2<CCTK_REAL> gf_Gamtz_rhs_(&layout,
                                       Gamtz_rhs);

  const GF3D2<CCTK_REAL> gf_Theta_rhs_(&layout,
                                       Theta_rhs);

  const GF3D2<CCTK_REAL> gf_alphaG_rhs_(&layout,
                                        alphaG_rhs);

  const GF3D2<CCTK_REAL> gf_betaGx_rhs_(&layout,
                                        betaGx_rhs);
  const GF3D2<CCTK_REAL> gf_betaGy_rhs_(&layout,
                                        betaGy_rhs);
  const GF3D2<CCTK_REAL> gf_betaGz_rhs_(&layout,
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
