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
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};

  const GF3D1<const CCTK_REAL> gf_chi_(cctkGH, indextype, nghostzones, chi);

  const GF3D1<const CCTK_REAL> gf_gammatxx_(cctkGH, indextype, nghostzones,
                                            gammatxx);
  const GF3D1<const CCTK_REAL> gf_gammatxy_(cctkGH, indextype, nghostzones,
                                            gammatxy);
  const GF3D1<const CCTK_REAL> gf_gammatxz_(cctkGH, indextype, nghostzones,
                                            gammatxz);
  const GF3D1<const CCTK_REAL> gf_gammatyy_(cctkGH, indextype, nghostzones,
                                            gammatyy);
  const GF3D1<const CCTK_REAL> gf_gammatyz_(cctkGH, indextype, nghostzones,
                                            gammatyz);
  const GF3D1<const CCTK_REAL> gf_gammatzz_(cctkGH, indextype, nghostzones,
                                            gammatzz);

  const GF3D1<const CCTK_REAL> gf_Kh_(cctkGH, indextype, nghostzones, Kh);

  const GF3D1<const CCTK_REAL> gf_Atxx_(cctkGH, indextype, nghostzones, Atxx);
  const GF3D1<const CCTK_REAL> gf_Atxy_(cctkGH, indextype, nghostzones, Atxy);
  const GF3D1<const CCTK_REAL> gf_Atxz_(cctkGH, indextype, nghostzones, Atxz);
  const GF3D1<const CCTK_REAL> gf_Atyy_(cctkGH, indextype, nghostzones, Atyy);
  const GF3D1<const CCTK_REAL> gf_Atyz_(cctkGH, indextype, nghostzones, Atyz);
  const GF3D1<const CCTK_REAL> gf_Atzz_(cctkGH, indextype, nghostzones, Atzz);

  const GF3D1<const CCTK_REAL> gf_Gamtx_(cctkGH, indextype, nghostzones, Gamtx);
  const GF3D1<const CCTK_REAL> gf_Gamty_(cctkGH, indextype, nghostzones, Gamty);
  const GF3D1<const CCTK_REAL> gf_Gamtz_(cctkGH, indextype, nghostzones, Gamtz);

  const GF3D1<const CCTK_REAL> gf_Theta_(cctkGH, indextype, nghostzones, Theta);

  const GF3D1<const CCTK_REAL> gf_alphaG_(cctkGH, indextype, nghostzones,
                                          alphaG);

  const GF3D1<const CCTK_REAL> gf_betaGx_(cctkGH, indextype, nghostzones,
                                          betaGx);
  const GF3D1<const CCTK_REAL> gf_betaGy_(cctkGH, indextype, nghostzones,
                                          betaGy);
  const GF3D1<const CCTK_REAL> gf_betaGz_(cctkGH, indextype, nghostzones,
                                          betaGz);

  const GF3D1<const CCTK_REAL> gf_eTtt_(cctkGH, indextype, nghostzones, eTtt);

  const GF3D1<const CCTK_REAL> gf_eTtx_(cctkGH, indextype, nghostzones, eTtx);
  const GF3D1<const CCTK_REAL> gf_eTty_(cctkGH, indextype, nghostzones, eTty);
  const GF3D1<const CCTK_REAL> gf_eTtz_(cctkGH, indextype, nghostzones, eTtz);

  const GF3D1<const CCTK_REAL> gf_eTxx_(cctkGH, indextype, nghostzones, eTxx);
  const GF3D1<const CCTK_REAL> gf_eTxy_(cctkGH, indextype, nghostzones, eTxy);
  const GF3D1<const CCTK_REAL> gf_eTxz_(cctkGH, indextype, nghostzones, eTxz);
  const GF3D1<const CCTK_REAL> gf_eTyy_(cctkGH, indextype, nghostzones, eTyy);
  const GF3D1<const CCTK_REAL> gf_eTyz_(cctkGH, indextype, nghostzones, eTyz);
  const GF3D1<const CCTK_REAL> gf_eTzz_(cctkGH, indextype, nghostzones, eTzz);

  const GF3D1<CCTK_REAL> gf_gxx_(cctkGH, indextype, nghostzones, gxx);
  const GF3D1<CCTK_REAL> gf_gxy_(cctkGH, indextype, nghostzones, gxy);
  const GF3D1<CCTK_REAL> gf_gxz_(cctkGH, indextype, nghostzones, gxz);
  const GF3D1<CCTK_REAL> gf_gyy_(cctkGH, indextype, nghostzones, gyy);
  const GF3D1<CCTK_REAL> gf_gyz_(cctkGH, indextype, nghostzones, gyz);
  const GF3D1<CCTK_REAL> gf_gzz_(cctkGH, indextype, nghostzones, gzz);

  const GF3D1<CCTK_REAL> gf_kxx_(cctkGH, indextype, nghostzones, kxx);
  const GF3D1<CCTK_REAL> gf_kxy_(cctkGH, indextype, nghostzones, kxy);
  const GF3D1<CCTK_REAL> gf_kxz_(cctkGH, indextype, nghostzones, kxz);
  const GF3D1<CCTK_REAL> gf_kyy_(cctkGH, indextype, nghostzones, kyy);
  const GF3D1<CCTK_REAL> gf_kyz_(cctkGH, indextype, nghostzones, kyz);
  const GF3D1<CCTK_REAL> gf_kzz_(cctkGH, indextype, nghostzones, kzz);

  const GF3D1<CCTK_REAL> gf_alp_(cctkGH, indextype, nghostzones, alp);

  const GF3D1<CCTK_REAL> gf_dtalp_(cctkGH, indextype, nghostzones, dtalp);

  const GF3D1<CCTK_REAL> gf_betax_(cctkGH, indextype, nghostzones, betax);
  const GF3D1<CCTK_REAL> gf_betay_(cctkGH, indextype, nghostzones, betay);
  const GF3D1<CCTK_REAL> gf_betaz_(cctkGH, indextype, nghostzones, betaz);

  const GF3D1<CCTK_REAL> gf_dtbetax_(cctkGH, indextype, nghostzones, dtbetax);
  const GF3D1<CCTK_REAL> gf_dtbetay_(cctkGH, indextype, nghostzones, dtbetay);
  const GF3D1<CCTK_REAL> gf_dtbetaz_(cctkGH, indextype, nghostzones, dtbetaz);

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
