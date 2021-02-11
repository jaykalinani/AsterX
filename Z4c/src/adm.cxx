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

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_chi_(&layout, chi);

  const GF3D2<const CCTK_REAL> gf_gammatxx_(&layout, gammatxx);
  const GF3D2<const CCTK_REAL> gf_gammatxy_(&layout, gammatxy);
  const GF3D2<const CCTK_REAL> gf_gammatxz_(&layout, gammatxz);
  const GF3D2<const CCTK_REAL> gf_gammatyy_(&layout, gammatyy);
  const GF3D2<const CCTK_REAL> gf_gammatyz_(&layout, gammatyz);
  const GF3D2<const CCTK_REAL> gf_gammatzz_(&layout, gammatzz);

  const GF3D2<const CCTK_REAL> gf_Kh_(&layout, Kh);

  const GF3D2<const CCTK_REAL> gf_Atxx_(&layout, Atxx);
  const GF3D2<const CCTK_REAL> gf_Atxy_(&layout, Atxy);
  const GF3D2<const CCTK_REAL> gf_Atxz_(&layout, Atxz);
  const GF3D2<const CCTK_REAL> gf_Atyy_(&layout, Atyy);
  const GF3D2<const CCTK_REAL> gf_Atyz_(&layout, Atyz);
  const GF3D2<const CCTK_REAL> gf_Atzz_(&layout, Atzz);

  const GF3D2<const CCTK_REAL> gf_Gamtx_(&layout, Gamtx);
  const GF3D2<const CCTK_REAL> gf_Gamty_(&layout, Gamty);
  const GF3D2<const CCTK_REAL> gf_Gamtz_(&layout, Gamtz);

  const GF3D2<const CCTK_REAL> gf_Theta_(&layout, Theta);

  const GF3D2<const CCTK_REAL> gf_alphaG_(&layout, alphaG);

  const GF3D2<const CCTK_REAL> gf_betaGx_(&layout, betaGx);
  const GF3D2<const CCTK_REAL> gf_betaGy_(&layout, betaGy);
  const GF3D2<const CCTK_REAL> gf_betaGz_(&layout, betaGz);

  const GF3D2<const CCTK_REAL> gf_eTtt_(&layout, eTtt);

  const GF3D2<const CCTK_REAL> gf_eTtx_(&layout, eTtx);
  const GF3D2<const CCTK_REAL> gf_eTty_(&layout, eTty);
  const GF3D2<const CCTK_REAL> gf_eTtz_(&layout, eTtz);

  const GF3D2<const CCTK_REAL> gf_eTxx_(&layout, eTxx);
  const GF3D2<const CCTK_REAL> gf_eTxy_(&layout, eTxy);
  const GF3D2<const CCTK_REAL> gf_eTxz_(&layout, eTxz);
  const GF3D2<const CCTK_REAL> gf_eTyy_(&layout, eTyy);
  const GF3D2<const CCTK_REAL> gf_eTyz_(&layout, eTyz);
  const GF3D2<const CCTK_REAL> gf_eTzz_(&layout, eTzz);

  const GF3D2<CCTK_REAL> gf_gxx_(&layout, gxx);
  const GF3D2<CCTK_REAL> gf_gxy_(&layout, gxy);
  const GF3D2<CCTK_REAL> gf_gxz_(&layout, gxz);
  const GF3D2<CCTK_REAL> gf_gyy_(&layout, gyy);
  const GF3D2<CCTK_REAL> gf_gyz_(&layout, gyz);
  const GF3D2<CCTK_REAL> gf_gzz_(&layout, gzz);

  const GF3D2<CCTK_REAL> gf_kxx_(&layout, kxx);
  const GF3D2<CCTK_REAL> gf_kxy_(&layout, kxy);
  const GF3D2<CCTK_REAL> gf_kxz_(&layout, kxz);
  const GF3D2<CCTK_REAL> gf_kyy_(&layout, kyy);
  const GF3D2<CCTK_REAL> gf_kyz_(&layout, kyz);
  const GF3D2<CCTK_REAL> gf_kzz_(&layout, kzz);

  const GF3D2<CCTK_REAL> gf_alp_(&layout, alp);

  const GF3D2<CCTK_REAL> gf_dtalp_(&layout, dtalp);

  const GF3D2<CCTK_REAL> gf_betax_(&layout, betax);
  const GF3D2<CCTK_REAL> gf_betay_(&layout, betay);
  const GF3D2<CCTK_REAL> gf_betaz_(&layout, betaz);

  const GF3D2<CCTK_REAL> gf_dtbetax_(&layout, dtbetax);
  const GF3D2<CCTK_REAL> gf_dtbetay_(&layout, dtbetay);
  const GF3D2<CCTK_REAL> gf_dtbetaz_(&layout, dtbetaz);

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
        gf_eTtt_,                                                   //
        gf_eTtx_, gf_eTty_, gf_eTtz_,                               //
        gf_eTxx_, gf_eTxy_, gf_eTxz_, gf_eTyy_, gf_eTyz_, gf_eTzz_, //
        p.I);

    // Store
    vars.g.store(gf_gxx_, gf_gxy_, gf_gxz_, gf_gyy_, gf_gyz_, gf_gzz_, p.I);
    vars.K.store(gf_kxx_, gf_kxy_, gf_kxz_, gf_kyy_, gf_kyz_, gf_kzz_, p.I);
    gf_alp_(p.I) = vars.alp;
    gf_dtalp_(p.I) = vars.dtalp;
    vars.beta.store(gf_betax_, gf_betay_, gf_betaz_, p.I);
    vars.dtbeta.store(gf_dtbetax_, gf_dtbetay_, gf_dtbetaz_, p.I);
  });
}

} // namespace Z4c
