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
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_chi1(layout1, chi);

  const GF3D2<const CCTK_REAL> gf_gammatxx1(layout1, gammatxx);
  const GF3D2<const CCTK_REAL> gf_gammatxy1(layout1, gammatxy);
  const GF3D2<const CCTK_REAL> gf_gammatxz1(layout1, gammatxz);
  const GF3D2<const CCTK_REAL> gf_gammatyy1(layout1, gammatyy);
  const GF3D2<const CCTK_REAL> gf_gammatyz1(layout1, gammatyz);
  const GF3D2<const CCTK_REAL> gf_gammatzz1(layout1, gammatzz);

  const GF3D2<const CCTK_REAL> gf_Kh1(layout1, Kh);

  const GF3D2<const CCTK_REAL> gf_Atxx1(layout1, Atxx);
  const GF3D2<const CCTK_REAL> gf_Atxy1(layout1, Atxy);
  const GF3D2<const CCTK_REAL> gf_Atxz1(layout1, Atxz);
  const GF3D2<const CCTK_REAL> gf_Atyy1(layout1, Atyy);
  const GF3D2<const CCTK_REAL> gf_Atyz1(layout1, Atyz);
  const GF3D2<const CCTK_REAL> gf_Atzz1(layout1, Atzz);

  const GF3D2<const CCTK_REAL> gf_Gamtx1(layout1, Gamtx);
  const GF3D2<const CCTK_REAL> gf_Gamty1(layout1, Gamty);
  const GF3D2<const CCTK_REAL> gf_Gamtz1(layout1, Gamtz);

  const GF3D2<const CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<const CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const GF3D2<const CCTK_REAL> gf_betaGx1(layout1, betaGx);
  const GF3D2<const CCTK_REAL> gf_betaGy1(layout1, betaGy);
  const GF3D2<const CCTK_REAL> gf_betaGz1(layout1, betaGz);

  const GF3D2<const CCTK_REAL> gf_eTtt1(layout1, eTtt);

  const GF3D2<const CCTK_REAL> gf_eTtx1(layout1, eTtx);
  const GF3D2<const CCTK_REAL> gf_eTty1(layout1, eTty);
  const GF3D2<const CCTK_REAL> gf_eTtz1(layout1, eTtz);

  const GF3D2<const CCTK_REAL> gf_eTxx1(layout1, eTxx);
  const GF3D2<const CCTK_REAL> gf_eTxy1(layout1, eTxy);
  const GF3D2<const CCTK_REAL> gf_eTxz1(layout1, eTxz);
  const GF3D2<const CCTK_REAL> gf_eTyy1(layout1, eTyy);
  const GF3D2<const CCTK_REAL> gf_eTyz1(layout1, eTyz);
  const GF3D2<const CCTK_REAL> gf_eTzz1(layout1, eTzz);

  const GF3D2<CCTK_REAL> gf_gxx1(layout1, gxx);
  const GF3D2<CCTK_REAL> gf_gxy1(layout1, gxy);
  const GF3D2<CCTK_REAL> gf_gxz1(layout1, gxz);
  const GF3D2<CCTK_REAL> gf_gyy1(layout1, gyy);
  const GF3D2<CCTK_REAL> gf_gyz1(layout1, gyz);
  const GF3D2<CCTK_REAL> gf_gzz1(layout1, gzz);

  const GF3D2<CCTK_REAL> gf_kxx1(layout1, kxx);
  const GF3D2<CCTK_REAL> gf_kxy1(layout1, kxy);
  const GF3D2<CCTK_REAL> gf_kxz1(layout1, kxz);
  const GF3D2<CCTK_REAL> gf_kyy1(layout1, kyy);
  const GF3D2<CCTK_REAL> gf_kyz1(layout1, kyz);
  const GF3D2<CCTK_REAL> gf_kzz1(layout1, kzz);

  const GF3D2<CCTK_REAL> gf_alp1(layout1, alp);

  const GF3D2<CCTK_REAL> gf_dtalp1(layout1, dtalp);

  const GF3D2<CCTK_REAL> gf_betax1(layout1, betax);
  const GF3D2<CCTK_REAL> gf_betay1(layout1, betay);
  const GF3D2<CCTK_REAL> gf_betaz1(layout1, betaz);

  const GF3D2<CCTK_REAL> gf_dtbetax1(layout1, dtbetax);
  const GF3D2<CCTK_REAL> gf_dtbetay1(layout1, dtbetay);
  const GF3D2<CCTK_REAL> gf_dtbetaz1(layout1, dtbetaz);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    // Load and calculate
    const z4c_vars_noderivs<CCTK_REAL> vars(
        kappa1, kappa2, f_mu_L, f_mu_S, eta, //
        gf_chi1,                             //
        gf_gammatxx1, gf_gammatxy1, gf_gammatxz1, gf_gammatyy1, gf_gammatyz1,
        gf_gammatzz1,                                               //
        gf_Kh1,                                                     //
        gf_Atxx1, gf_Atxy1, gf_Atxz1, gf_Atyy1, gf_Atyz1, gf_Atzz1, //
        gf_Gamtx1, gf_Gamty1, gf_Gamtz1,                            //
        gf_Theta1,                                                  //
        gf_alphaG1,                                                 //
        gf_betaGx1, gf_betaGy1, gf_betaGz1,                         //
        gf_eTtt1,                                                   //
        gf_eTtx1, gf_eTty1, gf_eTtz1,                               //
        gf_eTxx1, gf_eTxy1, gf_eTxz1, gf_eTyy1, gf_eTyz1, gf_eTzz1, //
        p.I);

    // Store
    vars.g.store(gf_gxx1, gf_gxy1, gf_gxz1, gf_gyy1, gf_gyz1, gf_gzz1, p.I);
    vars.K.store(gf_kxx1, gf_kxy1, gf_kxz1, gf_kyy1, gf_kyz1, gf_kzz1, p.I);
    gf_alp1(p.I) = vars.alp;
    gf_dtalp1(p.I) = vars.dtalp;
    vars.beta.store(gf_betax1, gf_betay1, gf_betaz1, p.I);
    vars.dtbeta.store(gf_dtbetax1, gf_dtbetay1, gf_dtbetaz1, p.I);
  });
}

} // namespace Z4c
