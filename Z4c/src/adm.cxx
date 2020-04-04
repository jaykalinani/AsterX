#include "tensor.hxx"
#include "z4c_vars.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_ADM(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ADM;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_chi_(cctkGH, chi);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxx_(cctkGH, gammatxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxy_(cctkGH, gammatxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxz_(cctkGH, gammatxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatyy_(cctkGH, gammatyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatyz_(cctkGH, gammatyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatzz_(cctkGH, gammatzz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Kh_(cctkGH, Kh);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxx_(cctkGH, Atxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxy_(cctkGH, Atxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxz_(cctkGH, Atxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atyy_(cctkGH, Atyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atyz_(cctkGH, Atyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atzz_(cctkGH, Atzz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamtx_(cctkGH, Gamtx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamty_(cctkGH, Gamty);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamtz_(cctkGH, Gamtz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Theta_(cctkGH, Theta);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_alphaG_(cctkGH, alphaG);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGx_(cctkGH, betaGx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGy_(cctkGH, betaGy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGz_(cctkGH, betaGz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_gxx_(cctkGH, gxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gxy_(cctkGH, gxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gxz_(cctkGH, gxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gyy_(cctkGH, gyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gyz_(cctkGH, gyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gzz_(cctkGH, gzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_kxx_(cctkGH, kxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_kxy_(cctkGH, kxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_kxz_(cctkGH, kxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_kyy_(cctkGH, kyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_kyz_(cctkGH, kyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_kzz_(cctkGH, kzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_alp_(cctkGH, alp);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtalp_(cctkGH, dtalp);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_betax_(cctkGH, betax);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betay_(cctkGH, betay);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaz_(cctkGH, betaz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtbetax_(cctkGH, dtbetax);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtbetay_(cctkGH, dtbetay);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_dtbetaz_(cctkGH, dtbetaz);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    // Load and calculate
    const z4c_vars_noderivs<CCTK_REAL> vars(
        kappa1, kappa2, f_mu_L, f_mu_S, eta, //
        gf_chi_,                             //
        gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_, gf_gammatyz_,
        gf_gammatzz_,                                               //
        gf_Kh_,                                                     //
        gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_, gf_Atzz_, //
        gf_Gamtx_, gf_Gamty_, gf_Gamtz_,                            //
        gf_Theta_,                                                  //
        gf_alphaG_,                                                 //
        gf_betaGx_, gf_betaGy_, gf_betaGz_,                         //
        p.I);

    // Store
    vars.g.store(gf_gxx_, gf_gxy_, gf_gxz_, gf_gyy_, gf_gyz_, gf_gzz_, p.I);
    vars.k.store(gf_kxx_, gf_kxy_, gf_kxz_, gf_kyy_, gf_kyz_, gf_kzz_, p.I);
    gf_alp_(p.I) = vars.alp;
    vars.beta.store(gf_betax_, gf_betay_, gf_betaz_, p.I);
  });

#warning "TODO: These are wrong!"
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtalp_(p.I) = 0; });

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtbetax_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtbetay_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_dtbetaz_(p.I) = 0; });
}

} // namespace Z4c
