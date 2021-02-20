#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_Initial1(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Initial1;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout1(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_gxx1(layout1, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy1(layout1, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz1(layout1, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy1(layout1, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz1(layout1, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz1(layout1, gzz);

  const GF3D2<const CCTK_REAL> gf_kxx1(layout1, kxx);
  const GF3D2<const CCTK_REAL> gf_kxy1(layout1, kxy);
  const GF3D2<const CCTK_REAL> gf_kxz1(layout1, kxz);
  const GF3D2<const CCTK_REAL> gf_kyy1(layout1, kyy);
  const GF3D2<const CCTK_REAL> gf_kyz1(layout1, kyz);
  const GF3D2<const CCTK_REAL> gf_kzz1(layout1, kzz);

  const GF3D2<const CCTK_REAL> gf_alp1(layout1, alp);

  const GF3D2<const CCTK_REAL> gf_betax1(layout1, betax);
  const GF3D2<const CCTK_REAL> gf_betay1(layout1, betay);
  const GF3D2<const CCTK_REAL> gf_betaz1(layout1, betaz);

  const GF3D2<CCTK_REAL> gf_chi1(layout1, chi);

  const GF3D2<CCTK_REAL> gf_gammatxx1(layout1, gammatxx);
  const GF3D2<CCTK_REAL> gf_gammatxy1(layout1, gammatxy);
  const GF3D2<CCTK_REAL> gf_gammatxz1(layout1, gammatxz);
  const GF3D2<CCTK_REAL> gf_gammatyy1(layout1, gammatyy);
  const GF3D2<CCTK_REAL> gf_gammatyz1(layout1, gammatyz);
  const GF3D2<CCTK_REAL> gf_gammatzz1(layout1, gammatzz);

  const GF3D2<CCTK_REAL> gf_Kh1(layout1, Kh);

  const GF3D2<CCTK_REAL> gf_Atxx1(layout1, Atxx);
  const GF3D2<CCTK_REAL> gf_Atxy1(layout1, Atxy);
  const GF3D2<CCTK_REAL> gf_Atxz1(layout1, Atxz);
  const GF3D2<CCTK_REAL> gf_Atyy1(layout1, Atyy);
  const GF3D2<CCTK_REAL> gf_Atyz1(layout1, Atyz);
  const GF3D2<CCTK_REAL> gf_Atzz1(layout1, Atzz);

  const GF3D2<CCTK_REAL> gf_Theta1(layout1, Theta);

  const GF3D2<CCTK_REAL> gf_alphaG1(layout1, alphaG);

  const GF3D2<CCTK_REAL> gf_betaGx1(layout1, betaGx);
  const GF3D2<CCTK_REAL> gf_betaGy1(layout1, betaGy);
  const GF3D2<CCTK_REAL> gf_betaGz1(layout1, betaGz);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load
    const mat3<CCTK_REAL, DN, DN> g(gf_gxx1, gf_gxy1, gf_gxz1, gf_gyy1, gf_gyz1,
                                    gf_gzz1, p.I);
    const mat3<CCTK_REAL, DN, DN> k(gf_kxx1, gf_kxy1, gf_kxz1, gf_kyy1, gf_kyz1,
                                    gf_kzz1, p.I);
    const CCTK_REAL alp = gf_alp1(p.I);
    const vec3<CCTK_REAL, UP> beta(gf_betax1, gf_betay1, gf_betaz1, p.I);

    // Calculate Z4c variables (all except Gammat)
    const CCTK_REAL detg = g.det();
    const mat3<CCTK_REAL, UP, UP> gu = g.inv(detg);

    const CCTK_REAL chi = 1 / cbrt(detg);

    const mat3<CCTK_REAL, DN, DN> gammat(
        [&](int a, int b) { return chi * g(a, b); });

    const CCTK_REAL K = sum2([&](int x, int y) { return gu(x, y) * k(x, y); });

    const CCTK_REAL Theta = 0;

    const CCTK_REAL Kh = K - 2 * Theta;

    const mat3<CCTK_REAL, DN, DN> At(
        [&](int a, int b) { return chi * (k(a, b) - K / 3 * g(a, b)); });

    const CCTK_REAL alphaG = alp;

    const vec3<CCTK_REAL, UP> betaG([&](int a) { return beta(a); });

    // Store
    gf_chi1(p.I) = chi;
    gammat.store(gf_gammatxx1, gf_gammatxy1, gf_gammatxz1, gf_gammatyy1,
                 gf_gammatyz1, gf_gammatzz1, p.I);
    gf_Kh1(p.I) = Kh;
    At.store(gf_Atxx1, gf_Atxy1, gf_Atxz1, gf_Atyy1, gf_Atyz1, gf_Atzz1, p.I);
    gf_Theta1(p.I) = Theta;
    gf_alphaG1(p.I) = alphaG;
    betaG.store(gf_betaGx1, gf_betaGy1, gf_betaGz1, p.I);
  });
}

} // namespace Z4c
