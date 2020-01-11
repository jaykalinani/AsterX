#warning "TODO"
#include <iostream>

#include "physics.hxx"
#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_RHS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_RHS;
  DECLARE_CCTK_PARAMETERS;

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

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_alphaG_(cctkGH, alphaG);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGx_(cctkGH, betaGx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGy_(cctkGH, betaGy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_betaGz_(cctkGH, betaGz);

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

    const CCTK_REAL alphaG = gf_alphaG_(p.I);

    const vec3<CCTK_REAL> betaG(gf_betaGx_, gf_betaGy_, gf_betaGz_, p);

    constexpr CCTK_REAL rho = 0;
    constexpr vec3<CCTK_REAL> S{0, 0, 0};
    constexpr CCTK_REAL traceT = 0;
    constexpr mat3<CCTK_REAL> T{0, 0, 0, 0, 0, 0};

    // Calculate RHS
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
    mat3<vec3<CCTK_REAL> > dg;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        for (int c = 0; c < 3; ++c)
          dg(a, b)(c) = dchi(c) * g(a, b) + chi * dgammat(a, b)(c);

    const mat3<mat3<CCTK_REAL> > ddgammat{
        deriv2(gf_gammatxx_, p.I, dx), deriv2(gf_gammatxy_, p.I, dx),
        deriv2(gf_gammatxz_, p.I, dx), deriv2(gf_gammatyy_, p.I, dx),
        deriv2(gf_gammatyz_, p.I, dx), deriv2(gf_gammatzz_, p.I, dx),
    };

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

    vec3<mat3<CCTK_REAL> > Gammatl;
    vec3<mat3<CCTK_REAL> > Gammat;
    vec3<CCTK_REAL> Gamtd;
    calc_gamma(gammat, gammatu, dgammat, Gammatl, Gammat, Gamtd);
    vec3<mat3<CCTK_REAL> > Gammal;
    vec3<mat3<CCTK_REAL> > Gamma;
    vec3<CCTK_REAL> Gam;
    calc_gamma(g, gu, dg, Gammal, Gamma, Gam);

    mat3<CCTK_REAL> DDchi;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b) {
        CCTK_REAL s = ddchi(a, b);
        for (int x = 0; x < 3; ++x)
          s -= Gammat(x)(a, b) * dchi(x);
        DDchi(a, b) = s;
      }

    const vec3<vec3<CCTK_REAL> > dGamt{
        deriv(gf_Gamtx_, p.I, dx),
        deriv(gf_Gamty_, p.I, dx),
        deriv(gf_Gamtz_, p.I, dx),
    };

    const mat3<CCTK_REAL> R =
        calc_ricci(chi, dchi, DDchi, gammat, gammatu, ddgammat, Gammatl, Gammat,
                   Gamtd, dGamt);

    CCTK_REAL Rsc = 0;
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        Rsc += gu(x, y) * R(x, y);

    const vec3<CCTK_REAL> dalphaG = deriv(gf_alphaG_, p.I, dx);
    const mat3<CCTK_REAL> ddalphaG = deriv2(gf_alphaG_, p.I, dx);

    const vec3<CCTK_REAL> dTheta = deriv(gf_Theta_, p.I, dx);

    mat3<CCTK_REAL> DDalphaG;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b) {
        CCTK_REAL s = ddalphaG(a, b);
        for (int x = 0; x < 3; ++x)
          s -= Gammat(x)(a, b) * dalphaG(x);
        DDalphaG(a, b) = s;
      }

    const vec3<vec3<CCTK_REAL> > dbetaG{
        deriv(gf_betaGx_, p.I, dx),
        deriv(gf_betaGy_, p.I, dx),
        deriv(gf_betaGz_, p.I, dx),
    };

    const vec3<mat3<CCTK_REAL> > ddbetaG{
        deriv2(gf_betaGx_, p.I, dx),
        deriv2(gf_betaGy_, p.I, dx),
        deriv2(gf_betaGz_, p.I, dx),
    };

    // (1)
    CCTK_REAL chi_rhs = 0;
    chi_rhs += alphaG * (Kh + 2 * Theta);
    for (int x = 0; x < 3; ++x)
      chi_rhs -= dbetaG(x)(x);
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        chi_rhs -= Gamma(x)(x, y) * betaG(y);
    chi_rhs *= chi * 2 / 3;

    // (2)
    mat3<CCTK_REAL> gammat_rhs;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b) {
        CCTK_REAL s = 0;
        s += -2 * alphaG * At(a, b);
        for (int x = 0; x < 3; ++x)
          s += betaG(x) * dgammat(a, b)(x);
        for (int x = 0; x < 3; ++x)
          s += 2 * (gammat(x, a) * dbetaG(x)(b) + gammat(x, b) * dbetaG(x)(a));
        for (int x = 0; x < 3; ++x)
          s -= gammat(a, b) * dbetaG(x)(x) * 2 / 3;
        gammat_rhs(a, b) = s;
      }

    // (3)
    CCTK_REAL Kh_rhs = 0;
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        Kh_rhs -= gu(x, y) * DDalphaG(x, y);
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        Kh_rhs += alphaG * (At(x, y) * Atu(x, y) + pow2(Kh + 2 * Theta) / 3);
    Kh_rhs += 4 * M_PI * alphaG * (traceT + rho);
    Kh_rhs += alphaG * kappa1 * (1 - kappa2) * Theta;
    for (int x = 0; x < 3; ++x)
      Kh_rhs += betaG(x) * dKh(x);

    // (4)
    mat3<CCTK_REAL> At_rhs1;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b) {
        CCTK_REAL s = 0;
        s -= DDalphaG(a, b) + alphaG * (R(a, b) - 8 * M_PI * T(a, b));
        At_rhs1(a, b) = s;
      }
    CCTK_REAL trace_At_rhs1 = 0;
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        trace_At_rhs1 += gammatu(x, y) * At_rhs1(x, y);
    mat3<CCTK_REAL> At_rhs;
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b) {
        CCTK_REAL s = 0;
        s += chi * (At_rhs1(a, b) - 1 / 3 * trace_At_rhs1 * gammat(a, b));
        s += alphaG * (Kh + 2 * Theta) * At(a, b);
        for (int x = 0; x < 3; ++x)
          for (int y = 0; y < 3; ++y)
            s -= alphaG * 2 * gammatu(x, y) * At(x, a) * At(y, b);
        for (int x = 0; x < 3; ++x)
          s += betaG(x) * dAt(a, b)(x);
        for (int x = 0; x < 3; ++x)
          s += 2 * (At(x, a) * dbetaG(x)(b) + At(x, b) * dbetaG(x)(a));
        for (int x = 0; x < 3; ++x)
          s -= At(a, b) * dbetaG(x)(x) * 2 / 3;
        At_rhs(a, b) = s;
      }

    // (5)
    vec3<CCTK_REAL> Gamt_rhs;
    for (int a = 0; a < 3; ++a) {
      CCTK_REAL s = 0;
      for (int x = 0; x < 3; ++x)
        s -= 2 * Atu(a, x) * dalphaG(x);
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
          s += 2 * alphaG * Gammat(a)(x, y) * Atu(x, y);
      for (int x = 0; x < 3; ++x)
        s -= 2 * alphaG * Atu(a, x) * dchi(x) / chi * 3 / 2;
      for (int x = 0; x < 3; ++x)
        s -= 2 * alphaG * gammatu(a, x) * (2 * dKh(x) + dTheta(x)) / 3;
      for (int x = 0; x < 3; ++x)
        s -= 2 * alphaG * 8 * M_PI * gammatu(a, x) * S(x);
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
          s += gammatu(x, y) * ddbetaG(a)(x, y);
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
          s += gammatu(a, x) * ddbetaG(y)(x, y) / 3;
      for (int x = 0; x < 3; ++x)
        s += betaG(x) * dGamt(a)(x);
      for (int x = 0; x < 3; ++x)
        s -= Gamtd(x) * dbetaG(a)(x);
      for (int x = 0; x < 3; ++x)
        s += Gamtd(a) * dbetaG(x)(x) * 2 / 3;
      s -= 2 * alphaG * kappa1 * (Gamt(a) - Gamtd(a));
      Gamt_rhs(a) = s;
    }

    // (6)
    CCTK_REAL Theta_rhs = 0;
    Theta_rhs += alphaG * Rsc / 2;
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        Theta_rhs -= alphaG * At(x, y) * Atu(x, y) / 2;
    Theta_rhs += alphaG * pow2(Kh + 2 * Theta) / 3;
    Theta_rhs -= alphaG * 8 * M_PI * rho;
    Theta_rhs -= alphaG * kappa1 * (2 + kappa2) * Theta;
    for (int x = 0; x < 3; ++x)
      Theta_rhs += betaG(x) * dTheta(x);

    const CCTK_REAL mu_L = f_mu_L / alphaG;

    CCTK_REAL alphaG_rhs = 0;
    alphaG_rhs -= pow2(alphaG) * mu_L * Kh;
    for (int x = 0; x < 3; ++x)
      alphaG_rhs += betaG(x) * dalphaG(x);

    const CCTK_REAL mu_S = f_mu_S / pow2(alphaG);

    vec3<CCTK_REAL> betaG_rhs;
    for (int a = 0; a < 3; ++a) {
      CCTK_REAL s = 0;
      s += pow2(alphaG) * mu_S * Gamt(a) - eta * betaG(a);
      for (int x = 0; x < 3; ++x)
        s += betaG(x) * dbetaG(a)(x);
      betaG_rhs(a) = s;
    }

    // Store
    gf_chi_rhs_(p.I) = chi_rhs;
    gammat_rhs.store(gf_gammatxx_rhs_, gf_gammatxy_rhs_, gf_gammatxz_rhs_,
                     gf_gammatyy_rhs_, gf_gammatyz_rhs_, gf_gammatzz_rhs_, p);
    gf_Kh_rhs_(p.I) = Kh_rhs;
    At_rhs.store(gf_Atxx_rhs_, gf_Atxy_rhs_, gf_Atxz_rhs_, gf_Atyy_rhs_,
                 gf_Atyz_rhs_, gf_Atzz_rhs_, p);
    Gamt_rhs.store(gf_Gamtx_rhs_, gf_Gamty_rhs_, gf_Gamtz_rhs_, p);
    gf_Theta_rhs_(p.I) = Theta_rhs;
    gf_alphaG_rhs_(p.I) = alphaG_rhs;
    betaG_rhs.store(gf_betaGx_rhs_, gf_betaGy_rhs_, gf_betaGz_rhs_, p);
  });
}

} // namespace Z4c
