#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHSBoundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_chi_rhs1(layout1, chi_rhs);

  const GF3D2<CCTK_REAL> gf_gammatxx_rhs1(layout1, gammatxx_rhs);
  const GF3D2<CCTK_REAL> gf_gammatxy_rhs1(layout1, gammatxy_rhs);
  const GF3D2<CCTK_REAL> gf_gammatxz_rhs1(layout1, gammatxz_rhs);
  const GF3D2<CCTK_REAL> gf_gammatyy_rhs1(layout1, gammatyy_rhs);
  const GF3D2<CCTK_REAL> gf_gammatyz_rhs1(layout1, gammatyz_rhs);
  const GF3D2<CCTK_REAL> gf_gammatzz_rhs1(layout1, gammatzz_rhs);

  const GF3D2<CCTK_REAL> gf_Kh_rhs1(layout1, Kh_rhs);

  const GF3D2<CCTK_REAL> gf_Atxx_rhs1(layout1, Atxx_rhs);
  const GF3D2<CCTK_REAL> gf_Atxy_rhs1(layout1, Atxy_rhs);
  const GF3D2<CCTK_REAL> gf_Atxz_rhs1(layout1, Atxz_rhs);
  const GF3D2<CCTK_REAL> gf_Atyy_rhs1(layout1, Atyy_rhs);
  const GF3D2<CCTK_REAL> gf_Atyz_rhs1(layout1, Atyz_rhs);
  const GF3D2<CCTK_REAL> gf_Atzz_rhs1(layout1, Atzz_rhs);

  const GF3D2<CCTK_REAL> gf_Gamtx_rhs1(layout1, Gamtx_rhs);
  const GF3D2<CCTK_REAL> gf_Gamty_rhs1(layout1, Gamty_rhs);
  const GF3D2<CCTK_REAL> gf_Gamtz_rhs1(layout1, Gamtz_rhs);

  const GF3D2<CCTK_REAL> gf_Theta_rhs1(layout1, Theta_rhs);

  const GF3D2<CCTK_REAL> gf_alphaG_rhs1(layout1, alphaG_rhs);

  const GF3D2<CCTK_REAL> gf_betaGx_rhs1(layout1, betaGx_rhs);
  const GF3D2<CCTK_REAL> gf_betaGy_rhs1(layout1, betaGy_rhs);
  const GF3D2<CCTK_REAL> gf_betaGz_rhs1(layout1, betaGz_rhs);

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_chi_rhs1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatxx_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatxy_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatxz_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatyy_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatyz_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_gammatzz_rhs1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Kh_rhs1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxx_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxy_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atxz_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atyy_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atyz_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_Atzz_rhs1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_Gamtx_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_Gamty_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_Gamtz_rhs1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_Theta_rhs1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_alphaG_rhs1(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_betaGx_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_betaGy_rhs1(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { gf_betaGz_rhs1(p.I) = 0; });
}

} // namespace Z4c
