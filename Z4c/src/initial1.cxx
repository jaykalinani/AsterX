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

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gxx_(cctkGH, gxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gxy_(cctkGH, gxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gxz_(cctkGH, gxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gyy_(cctkGH, gyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gyz_(cctkGH, gyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gzz_(cctkGH, gzz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_kxx_(cctkGH, kxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_kxy_(cctkGH, kxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_kxz_(cctkGH, kxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_kyy_(cctkGH, kyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_kyz_(cctkGH, kyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_kzz_(cctkGH, kzz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_alp_(cctkGH, alp);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betax_(cctkGH, betax);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betay_(cctkGH, betay);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaz_(cctkGH, betaz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_chi_(cctkGH, chi);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxx_(cctkGH, gammatxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxy_(cctkGH, gammatxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatxz_(cctkGH, gammatxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyy_(cctkGH, gammatyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatyz_(cctkGH, gammatyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_gammatzz_(cctkGH, gammatzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Kh_(cctkGH, Kh);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxx_(cctkGH, Atxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxy_(cctkGH, Atxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atxz_(cctkGH, Atxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyy_(cctkGH, Atyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atyz_(cctkGH, Atyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_Atzz_(cctkGH, Atzz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_Theta_(cctkGH, Theta);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_alphaG_(cctkGH, alphaG);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGx_(cctkGH, betaGx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGy_(cctkGH, betaGy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_betaGz_(cctkGH, betaGz);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load
    const mat3<CCTK_REAL> g(gf_gxx_, gf_gxy_, gf_gxz_, gf_gyy_, gf_gyz_,
                            gf_gzz_, p.I);
    const mat3<CCTK_REAL> k(gf_kxx_, gf_kxy_, gf_kxz_, gf_kyy_, gf_kyz_,
                            gf_kzz_, p.I);
    const CCTK_REAL alp = gf_alp_(p.I);
    const vec3<CCTK_REAL> beta(gf_betax_, gf_betay_, gf_betaz_, p.I);

    // Calculate Z4c variables (all except Gammat)
    const CCTK_REAL detg = g.det();
    const mat3<CCTK_REAL> gu = g.inv(detg);

    const CCTK_REAL chi = 1 / cbrt(detg);

    const mat3<CCTK_REAL> gammat([&](int a, int b) { return chi * g(a, b); });

    const CCTK_REAL K = sum2([&](int x, int y) { return gu(x, y) * k(x, y); });

    const CCTK_REAL Theta = 0;

    const CCTK_REAL Kh = K - 2 * Theta;

    const mat3<CCTK_REAL> At(
        [&](int a, int b) { return chi * (k(a, b) - K / 3 * g(a, b)); });

    const CCTK_REAL alphaG = alp;

    const vec3<CCTK_REAL> betaG([&](int a) { return beta(a); });

    // Store
    gf_chi_(p.I) = chi;
    gammat.store(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
                 gf_gammatyz_, gf_gammatzz_, p);
    gf_Kh_(p.I) = Kh;
    At.store(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_, gf_Atzz_, p);
    gf_Theta_(p.I) = Theta;
    gf_alphaG_(p.I) = alphaG;
    betaG.store(gf_betaGx_, gf_betaGy_, gf_betaGz_, p);
  });
}

} // namespace Z4c
