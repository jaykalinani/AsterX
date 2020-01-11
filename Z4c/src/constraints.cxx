#include "physics.hxx"
#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Constraints;

  const vec3<CCTK_REAL> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

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

  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCx_(cctkGH, ZtCx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCy_(cctkGH, ZtCy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCz_(cctkGH, ZtCz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_HC_(cctkGH, HC);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCx_(cctkGH, MtCx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCy_(cctkGH, MtCy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCz_(cctkGH, MtCz);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load
    const CCTK_REAL chi = gf_chi_(p.I);
    const mat3<CCTK_REAL> gammat(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_,
                                 gf_gammatyy_, gf_gammatyz_, gf_gammatzz_, p);
    const CCTK_REAL Kh = gf_Kh_(p.I);
    const mat3<CCTK_REAL> At(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_, gf_Atyz_,
                             gf_Atzz_, p);
    const vec3<CCTK_REAL> Gamt(gf_Gamtx_, gf_Gamty_, gf_Gamtz_, p);
    const CCTK_REAL Theta = gf_Theta_(p.I);

    constexpr CCTK_REAL rho = 0;
    constexpr vec3<CCTK_REAL> S{0, 0, 0};

    // Calculate constraints
    // See arXiv:1212.2901 [gr-qc].

    const CCTK_REAL chi1 = 1 / chi;

    const vec3<CCTK_REAL> dchi = deriv(gf_chi_, p.I, dx);
    const mat3<CCTK_REAL> ddchi = deriv2(gf_chi_, p.I, dx);

    mat3<CCTK_REAL> g;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        g(a, b) = chi * gammat(a, b);

    const mat3<CCTK_REAL> gammatu = gammat.inv(1);
    mat3<CCTK_REAL> gu;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        gu(a, b) = chi1 * gammatu(a, b);

    const mat3<vec3<CCTK_REAL> > dgammat{
        deriv(gf_gammatxx_, p.I, dx), deriv(gf_gammatxy_, p.I, dx),
        deriv(gf_gammatxz_, p.I, dx), deriv(gf_gammatyy_, p.I, dx),
        deriv(gf_gammatyz_, p.I, dx), deriv(gf_gammatzz_, p.I, dx),
    };

    const mat3<mat3<CCTK_REAL> > ddgammat{
        deriv2(gf_gammatxx_, p.I, dx), deriv2(gf_gammatxy_, p.I, dx),
        deriv2(gf_gammatxz_, p.I, dx), deriv2(gf_gammatyy_, p.I, dx),
        deriv2(gf_gammatyz_, p.I, dx), deriv2(gf_gammatzz_, p.I, dx),
    };

    const mat3<vec3<CCTK_REAL> > dgammatu = calc_dgu(gammatu, dgammat);

    const vec3<CCTK_REAL> dKh = deriv(gf_Kh_, p.I, dx);

    mat3<CCTK_REAL> Atu;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b) {
        CCTK_REAL s = 0;
        for (int x = 0; x < 3; ++x)
          for (int y = 0; y < 3; ++y)
            s += gammatu(a, x) * gammatu(b, y) * At(x, y);
        Atu(a, b) = s;
      }

    const mat3<vec3<CCTK_REAL> > dAt{
        deriv(gf_Atxx_, p.I, dx), deriv(gf_Atxy_, p.I, dx),
        deriv(gf_Atxz_, p.I, dx), deriv(gf_Atyy_, p.I, dx),
        deriv(gf_Atyz_, p.I, dx), deriv(gf_Atzz_, p.I, dx),
    };

    const mat3<vec3<CCTK_REAL> > dAtu = calc_dAu(gammatu, dgammatu, At, dAt);

    const vec3<vec3<CCTK_REAL> > dGamt{
        deriv(gf_Gamtx_, p.I, dx),
        deriv(gf_Gamty_, p.I, dx),
        deriv(gf_Gamtz_, p.I, dx),
    };

    vec3<mat3<CCTK_REAL> > Gammatdl;
    vec3<mat3<CCTK_REAL> > Gammatd;
    vec3<CCTK_REAL> Gamtd;
    calc_gamma(gammat, gammatu, dgammat, Gammatdl, Gammatd, Gamtd);

    mat3<CCTK_REAL> DDchi;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b) {
        CCTK_REAL s = ddchi(a, b);
        for (int x = 0; x < 3; ++x)
          s -= Gammatd(x)(a, b) * dchi(x);
        DDchi(a, b) = s;
      }

    const vec3<CCTK_REAL> dTheta = deriv(gf_Theta_, p.I, dx);

    // (13)
    vec3<CCTK_REAL> ZtC;
    for (int a = 0; a < 3; ++a)
      ZtC(a) = (Gamt(a) - Gamtd(a)) / 2;

    const mat3<CCTK_REAL> R =
        calc_ricci(chi, dchi, DDchi, gammat, gammatu, ddgammat, Gammatdl,
                   Gammatd, Gamtd, dGamt);

    CCTK_REAL Rsc = 0;
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        Rsc += gu(x, y) * R(x, y);

    // (14)
    CCTK_REAL HC = 0;
    HC += Rsc;
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        HC += At(x, y) * Atu(x, y);
    HC -= pow2(Kh + 2 * Theta) * 2 / 3;
    HC -= 16 * M_PI * rho;

    // (15)
    vec3<CCTK_REAL> MtC;
    for (int a = 0; a < 3; ++a) {
      CCTK_REAL s = 0;
      for (int x = 0; x < 3; ++x)
        s += dAtu(a, x)(x);
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
          s += Gammatd(a)(x, y) * Atu(x, y);
      for (int x = 0; x < 3; ++x)
        s -= gammatu(a, x) * (dKh(x) + 2 * dTheta(x)) * 2 / 3;
      for (int x = 0; x < 3; ++x)
        s -= Atu(a, x) * dchi(x) / chi * 2 / 3;
      for (int x = 0; x < 3; ++x)
        s -= 8 * M_PI * gammatu(a, x) * S(x);
      MtC(a) = s;
    }

    // Store
    ZtC.store(gf_ZtCx_, gf_ZtCy_, gf_ZtCz_, p);
    gf_HC_(p.I) = HC;
    MtC.store(gf_MtCx_, gf_MtCy_, gf_MtCz_, p);
  });
}

} // namespace Z4c
