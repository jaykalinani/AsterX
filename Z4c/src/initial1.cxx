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
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<const CCTK_REAL> gf_gxx_(&layout, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy_(&layout, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz_(&layout, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy_(&layout, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz_(&layout, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz_(&layout, gzz);

  const GF3D2<const CCTK_REAL> gf_kxx_(&layout, kxx);
  const GF3D2<const CCTK_REAL> gf_kxy_(&layout, kxy);
  const GF3D2<const CCTK_REAL> gf_kxz_(&layout, kxz);
  const GF3D2<const CCTK_REAL> gf_kyy_(&layout, kyy);
  const GF3D2<const CCTK_REAL> gf_kyz_(&layout, kyz);
  const GF3D2<const CCTK_REAL> gf_kzz_(&layout, kzz);

  const GF3D2<const CCTK_REAL> gf_alp_(&layout, alp);

  const GF3D2<const CCTK_REAL> gf_betax_(&layout, betax);
  const GF3D2<const CCTK_REAL> gf_betay_(&layout, betay);
  const GF3D2<const CCTK_REAL> gf_betaz_(&layout, betaz);

  const GF3D2<CCTK_REAL> gf_chi_(&layout, chi);

  const GF3D2<CCTK_REAL> gf_gammatxx_(&layout, gammatxx);
  const GF3D2<CCTK_REAL> gf_gammatxy_(&layout, gammatxy);
  const GF3D2<CCTK_REAL> gf_gammatxz_(&layout, gammatxz);
  const GF3D2<CCTK_REAL> gf_gammatyy_(&layout, gammatyy);
  const GF3D2<CCTK_REAL> gf_gammatyz_(&layout, gammatyz);
  const GF3D2<CCTK_REAL> gf_gammatzz_(&layout, gammatzz);

  const GF3D2<CCTK_REAL> gf_Kh_(&layout, Kh);

  const GF3D2<CCTK_REAL> gf_Atxx_(&layout, Atxx);
  const GF3D2<CCTK_REAL> gf_Atxy_(&layout, Atxy);
  const GF3D2<CCTK_REAL> gf_Atxz_(&layout, Atxz);
  const GF3D2<CCTK_REAL> gf_Atyy_(&layout, Atyy);
  const GF3D2<CCTK_REAL> gf_Atyz_(&layout, Atyz);
  const GF3D2<CCTK_REAL> gf_Atzz_(&layout, Atzz);

  const GF3D2<CCTK_REAL> gf_Theta_(&layout, Theta);

  const GF3D2<CCTK_REAL> gf_alphaG_(&layout, alphaG);

  const GF3D2<CCTK_REAL> gf_betaGx_(&layout, betaGx);
  const GF3D2<CCTK_REAL> gf_betaGy_(&layout, betaGy);
  const GF3D2<CCTK_REAL> gf_betaGz_(&layout, betaGz);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load
    const mat3<CCTK_REAL, DN, DN> g(gf_gxx_, gf_gxy_, gf_gxz_, gf_gyy_, gf_gyz_,
                                    gf_gzz_, p.I);
    const mat3<CCTK_REAL, DN, DN> k(gf_kxx_, gf_kxy_, gf_kxz_, gf_kyy_, gf_kyz_,
                                    gf_kzz_, p.I);
    const CCTK_REAL alp = gf_alp_(p.I);
    const vec3<CCTK_REAL, UP> beta(gf_betax_, gf_betay_, gf_betaz_, p.I);

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
    gf_chi_(p.I) = chi;
    gammat.store(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
                 gf_gammatyz_, gf_gammatzz_, p.I);
    gf_Kh_(p.I) = Kh;
    At.store(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_, gf_Atzz_, p.I);
    gf_Theta_(p.I) = Theta;
    gf_alphaG_(p.I) = alphaG;
    betaG.store(gf_betaGx_, gf_betaGy_, gf_betaGz_, p.I);
  });
}

} // namespace Z4c
