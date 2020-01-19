#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_RHSBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHSBoundaries;

  const GF3D<CCTK_REAL, 0, 0, 0> gf_chi_rhs_(cctkGH, chi_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxx_rhs_(cctkGH, gammatxx_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxy_rhs_(cctkGH, gammatxy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxz_rhs_(cctkGH, gammatxz_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyy_rhs_(cctkGH, gammatyy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyz_rhs_(cctkGH, gammatyz_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatzz_rhs_(cctkGH, gammatzz_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Kh_rhs_(cctkGH, Kh_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxx_rhs_(cctkGH, Atxx_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxy_rhs_(cctkGH, Atxy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxz_rhs_(cctkGH, Atxz_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyy_rhs_(cctkGH, Atyy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyz_rhs_(cctkGH, Atyz_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atzz_rhs_(cctkGH, Atzz_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamtx_rhs_(cctkGH, Gamtx_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamty_rhs_(cctkGH, Gamty_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Gamtz_rhs_(cctkGH, Gamtz_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Theta_rhs_(cctkGH, Theta_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_alphaG_rhs_(cctkGH, alphaG_rhs);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGx_rhs_(cctkGH, betaGx_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGy_rhs_(cctkGH, betaGy_rhs);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGz_rhs_(cctkGH, betaGz_rhs);

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
