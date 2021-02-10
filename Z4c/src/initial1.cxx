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
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};

  const GF3D1<const CCTK_REAL> gf_gxx_(cctkGH, indextype, nghostzones, gxx);
  const GF3D1<const CCTK_REAL> gf_gxy_(cctkGH, indextype, nghostzones, gxy);
  const GF3D1<const CCTK_REAL> gf_gxz_(cctkGH, indextype, nghostzones, gxz);
  const GF3D1<const CCTK_REAL> gf_gyy_(cctkGH, indextype, nghostzones, gyy);
  const GF3D1<const CCTK_REAL> gf_gyz_(cctkGH, indextype, nghostzones, gyz);
  const GF3D1<const CCTK_REAL> gf_gzz_(cctkGH, indextype, nghostzones, gzz);

  const GF3D1<const CCTK_REAL> gf_kxx_(cctkGH, indextype, nghostzones, kxx);
  const GF3D1<const CCTK_REAL> gf_kxy_(cctkGH, indextype, nghostzones, kxy);
  const GF3D1<const CCTK_REAL> gf_kxz_(cctkGH, indextype, nghostzones, kxz);
  const GF3D1<const CCTK_REAL> gf_kyy_(cctkGH, indextype, nghostzones, kyy);
  const GF3D1<const CCTK_REAL> gf_kyz_(cctkGH, indextype, nghostzones, kyz);
  const GF3D1<const CCTK_REAL> gf_kzz_(cctkGH, indextype, nghostzones, kzz);

  const GF3D1<const CCTK_REAL> gf_alp_(cctkGH, indextype, nghostzones, alp);

  const GF3D1<const CCTK_REAL> gf_betax_(cctkGH, indextype, nghostzones, betax);
  const GF3D1<const CCTK_REAL> gf_betay_(cctkGH, indextype, nghostzones, betay);
  const GF3D1<const CCTK_REAL> gf_betaz_(cctkGH, indextype, nghostzones, betaz);

  const GF3D1<CCTK_REAL> gf_chi_(cctkGH, indextype, nghostzones, chi);

  const GF3D1<CCTK_REAL> gf_gammatxx_(cctkGH, indextype, nghostzones, gammatxx);
  const GF3D1<CCTK_REAL> gf_gammatxy_(cctkGH, indextype, nghostzones, gammatxy);
  const GF3D1<CCTK_REAL> gf_gammatxz_(cctkGH, indextype, nghostzones, gammatxz);
  const GF3D1<CCTK_REAL> gf_gammatyy_(cctkGH, indextype, nghostzones, gammatyy);
  const GF3D1<CCTK_REAL> gf_gammatyz_(cctkGH, indextype, nghostzones, gammatyz);
  const GF3D1<CCTK_REAL> gf_gammatzz_(cctkGH, indextype, nghostzones, gammatzz);

  const GF3D1<CCTK_REAL> gf_Kh_(cctkGH, indextype, nghostzones, Kh);

  const GF3D1<CCTK_REAL> gf_Atxx_(cctkGH, indextype, nghostzones, Atxx);
  const GF3D1<CCTK_REAL> gf_Atxy_(cctkGH, indextype, nghostzones, Atxy);
  const GF3D1<CCTK_REAL> gf_Atxz_(cctkGH, indextype, nghostzones, Atxz);
  const GF3D1<CCTK_REAL> gf_Atyy_(cctkGH, indextype, nghostzones, Atyy);
  const GF3D1<CCTK_REAL> gf_Atyz_(cctkGH, indextype, nghostzones, Atyz);
  const GF3D1<CCTK_REAL> gf_Atzz_(cctkGH, indextype, nghostzones, Atzz);

  const GF3D1<CCTK_REAL> gf_Theta_(cctkGH, indextype, nghostzones, Theta);

  const GF3D1<CCTK_REAL> gf_alphaG_(cctkGH, indextype, nghostzones, alphaG);

  const GF3D1<CCTK_REAL> gf_betaGx_(cctkGH, indextype, nghostzones, betaGx);
  const GF3D1<CCTK_REAL> gf_betaGy_(cctkGH, indextype, nghostzones, betaGy);
  const GF3D1<CCTK_REAL> gf_betaGz_(cctkGH, indextype, nghostzones, betaGz);

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
