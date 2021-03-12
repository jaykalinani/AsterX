#ifndef Z4C_HXX
#define Z4C_HXX

#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"

#include <cmath>
#include <iostream>

namespace Z4c {

// See arXiv:1212.2901 [gr-qc]
// Note: A_(ij) there means 1/2 (A_ij + A_ji)

template <typename T> struct z4c_vars_noderivs {

  // Parameters
  const T kappa1;
  const T kappa2;
  const T f_mu_L;
  const T f_mu_S;
  const T eta;

  // Tensor densities: TD = (sqrt detg)^W T
  //                   W[sqrt detg] = +1

  // Z4c variables
  const T chi;                  // W = -2/3
  const mat3<T, DN, DN> gammat; // W = -2/3
  const T Kh;                   // W = 0
  const mat3<T, DN, DN> At;     // W = -2/3
  const vec3<T, UP> Gamt;       // W = +2/3
  const T Theta;                // W = 0
  const T alphaG;               // W = 0
  const vec3<T, UP> betaG;      // W = 0

  // T_munu variables
  const T eTtt;
  const vec3<T, DN> eTti;
  const mat3<T, DN, DN> eTij;

  // Hydro variables
  const T rho;
  const vec3<T, DN> Si;
  const mat3<T, DN, DN> Sij;

  // ADM variables
  const mat3<T, DN, DN> g;  // W = 0
  const mat3<T, DN, DN> K;  // W = 0
  const T alp;              // W = 0
  const T dtalp;            // W = 0
  const vec3<T, UP> beta;   // W = 0
  const vec3<T, UP> dtbeta; // W = 0

  friend CCTK_ATTRIBUTE_NOINLINE ostream &
  operator<<(ostream &os, const z4c_vars_noderivs &vars) {
    return os << "z4c_vars_noderivs{"            //
              << "kappa1:" << vars.kappa1 << "," //
              << "kappa2:" << vars.kappa2 << "," //
              << "f_mu_L:" << vars.f_mu_L << "," //
              << "f_mu_S:" << vars.f_mu_S << "," //
              << "eta:" << vars.eta << ","       //
              << "chi:" << vars.chi << ","       //
              << "gammat:" << vars.gammat << "," //
              << "Kh:" << vars.Kh << ","         //
              << "At:" << vars.At << ","         //
              << "Gamt:" << vars.Gamt << ","     //
              << "Theta:" << vars.Theta << ","   //
              << "alphaG:" << vars.alphaG << "," //
              << "betaG:" << vars.betaG << ","   //
              << "eTtt:" << vars.eTtt << ","     //
              << "eTti:" << vars.eTti << ","     //
              << "eTij:" << vars.eTij << ","     //
              << "rho:" << vars.rho << ","       //
              << "Si:" << vars.Si << ","         //
              << "Sij:" << vars.Sij << ","       //
              << "g:" << vars.g << ","           //
              << "K:" << vars.K << ","           //
              << "alp:" << vars.alp << ","       //
              << "dtalp:" << vars.dtalp << ","   //
              << "dtbeta:" << vars.dtbeta << "," //
              << "}";
  }

  Z4C_INLINE Z4C_GPU z4c_vars_noderivs(
      const T &kappa1, const T &kappa2, const T &f_mu_L, const T &f_mu_S,
      const T &eta,
      //
      const T &chi, const mat3<T, DN, DN> &gammat, const T &Kh,
      const mat3<T, DN, DN> &At, const vec3<T, UP> &Gamt, const T &Theta,
      const T &alphaG, const vec3<T, UP> &betaG,
      //
      const T &eTtt, const vec3<T, DN> &eTti, const mat3<T, DN, DN> &eTij)
      : kappa1(kappa1), kappa2(kappa2), f_mu_L(f_mu_L), f_mu_S(f_mu_S),
        eta(eta),
        //
        chi(chi), gammat(gammat), Kh(Kh), At(At), Gamt(Gamt), Theta(Theta),
        alphaG(alphaG), betaG(betaG),
        //
        eTtt(eTtt), eTti(eTti), eTij(eTij),
        // Hydro variables
        // rho = n^a n^b T_ab
        rho(1 / pow2(alphaG) *
            (eTtt                                                             //
             - 2 * sum1([&] Z4C_INLINE(int x) { return betaG(x) * eTti(x); }) //
             + sum2([&] Z4C_INLINE(int x, int y) {
                 return betaG(x) * betaG(y) * eTij(x, y);
               }))),
        // S_i = -p_i^a n^b T_ab
        Si([&] Z4C_INLINE(int a) {
          return -1 / alphaG *
                 (eTti(a) //
                  - sum1([&] Z4C_INLINE(int x) {
                      return betaG(x) * eTij(a, x);
                    }));
        }), //
        // S_ij = p_i^a p_j^b T_ab
        Sij(eTij),
        // ADM variables
        g([&] Z4C_INLINE(int a, int b) { return 1 / chi * gammat(a, b); }), //
        K([&] Z4C_INLINE(int a, int b) {
          return 1 / chi * (At(a, b) + (Kh + 2 * Theta) / 3 * gammat(a, b));
        }),          //
        alp(alphaG), //
        // (11)
        dtalp([&] Z4C_INLINE {
          const T mu_L = f_mu_L / alphaG;
          return -pow2(alphaG) * mu_L * Kh;
        }()),
        beta([&] Z4C_INLINE(int a) { return betaG(a); }),
        // (12)
        dtbeta([&] Z4C_INLINE(int a) {
          const T mu_S = f_mu_S / pow2(alphaG);
          return pow2(alphaG) * mu_S * Gamt(a) //
                 - eta * betaG(a);
        })
  //
  {}

  Z4C_INLINE Z4C_GPU z4c_vars_noderivs(
      const T &kappa1, const T &kappa2, const T &f_mu_L, const T &f_mu_S,
      const T &eta,
      //
      const GF3D2<const T> &gf_chi_,
      //
      const GF3D2<const T> &gf_gammatxx_, const GF3D2<const T> &gf_gammatxy_,
      const GF3D2<const T> &gf_gammatxz_, const GF3D2<const T> &gf_gammatyy_,
      const GF3D2<const T> &gf_gammatyz_, const GF3D2<const T> &gf_gammatzz_,
      //
      const GF3D2<const T> &gf_Kh_,
      //
      const GF3D2<const T> &gf_Atxx_, const GF3D2<const T> &gf_Atxy_,
      const GF3D2<const T> &gf_Atxz_, const GF3D2<const T> &gf_Atyy_,
      const GF3D2<const T> &gf_Atyz_, const GF3D2<const T> &gf_Atzz_,
      //
      const GF3D2<const T> &gf_Gamtx_, const GF3D2<const T> &gf_Gamty_,
      const GF3D2<const T> &gf_Gamtz_,
      //
      const GF3D2<const T> &gf_Theta_,
      //
      const GF3D2<const T> &gf_alphaG_,
      //
      const GF3D2<const T> &gf_betaGx_, const GF3D2<const T> &gf_betaGy_,
      const GF3D2<const T> &gf_betaGz_,
      //
      const GF3D2<const T> &gf_eTtt_,
      //
      const GF3D2<const T> &gf_eTtx_, const GF3D2<const T> &gf_eTty_,
      const GF3D2<const T> &gf_eTtz_,
      //
      const GF3D2<const T> &gf_eTxx_, const GF3D2<const T> &gf_eTxy_,
      const GF3D2<const T> &gf_eTxz_, const GF3D2<const T> &gf_eTyy_,
      const GF3D2<const T> &gf_eTyz_, const GF3D2<const T> &gf_eTzz_,
      //
      const vect<int, 3> &I)
      : z4c_vars_noderivs<T>(kappa1, kappa2, f_mu_L, f_mu_S, eta,
                             //
                             gf_chi_(I),
                             //
                             mat3<T, DN, DN>(gf_gammatxx_, gf_gammatxy_,
                                             gf_gammatxz_, gf_gammatyy_,
                                             gf_gammatyz_, gf_gammatzz_, I),
                             //
                             gf_Kh_(I),
                             //
                             mat3<T, DN, DN>(gf_Atxx_, gf_Atxy_, gf_Atxz_,
                                             gf_Atyy_, gf_Atyz_, gf_Atzz_, I),
                             //
                             vec3<T, UP>(gf_Gamtx_, gf_Gamty_, gf_Gamtz_, I),
                             //
                             gf_Theta_(I),
                             //
                             gf_alphaG_(I),
                             //
                             vec3<T, UP>(gf_betaGx_, gf_betaGy_, gf_betaGz_, I),
                             //
                             gf_eTtt_(I),
                             //
                             vec3<T, UP>(gf_eTtx_, gf_eTty_, gf_eTtz_, I),
                             //
                             mat3<T, DN, DN>(gf_eTxx_, gf_eTxy_, gf_eTxz_,
                                             gf_eTyy_, gf_eTyz_, gf_eTzz_, I))
  //
  {}
};

template <typename T> struct z4c_vars : z4c_vars_noderivs<T> {

  // C++ is tedious:

  // Parameters
  using z4c_vars_noderivs<T>::kappa1;
  using z4c_vars_noderivs<T>::kappa2;
  using z4c_vars_noderivs<T>::f_mu_L;
  using z4c_vars_noderivs<T>::f_mu_S;
  using z4c_vars_noderivs<T>::eta;

  // Z4c variables
  using z4c_vars_noderivs<T>::chi;
  using z4c_vars_noderivs<T>::gammat;
  using z4c_vars_noderivs<T>::Kh;
  using z4c_vars_noderivs<T>::At;
  using z4c_vars_noderivs<T>::Gamt;
  using z4c_vars_noderivs<T>::Theta;
  using z4c_vars_noderivs<T>::alphaG;
  using z4c_vars_noderivs<T>::betaG;

  // T_munu variables
  using z4c_vars_noderivs<T>::eTtt;
  using z4c_vars_noderivs<T>::eTti;
  using z4c_vars_noderivs<T>::eTij;

  // Hydro variables
  using z4c_vars_noderivs<T>::rho;
  using z4c_vars_noderivs<T>::Si;
  using z4c_vars_noderivs<T>::Sij;

  // ADM variables
  using z4c_vars_noderivs<T>::g;
  using z4c_vars_noderivs<T>::K;
  using z4c_vars_noderivs<T>::alp;
  using z4c_vars_noderivs<T>::beta;
  using z4c_vars_noderivs<T>::dtalp;
  using z4c_vars_noderivs<T>::dtbeta;

  // Derivatives of Z4c variables
  const vec3<T, DN> dchi;
  const mat3<T, DN, DN> ddchi;
  const mat3<vec3<T, DN>, DN, DN> dgammat;
  const mat3<mat3<T, DN, DN>, DN, DN> ddgammat;
  const vec3<T, DN> dKh;
  const mat3<vec3<T, DN>, DN, DN> dAt;
  const vec3<vec3<T, DN>, UP> dGamt;
  const vec3<T, DN> dTheta;
  const vec3<T, DN> dalphaG;
  const mat3<T, DN, DN> ddalphaG;
  const vec3<vec3<T, DN>, UP> dbetaG;
  const vec3<mat3<T, DN, DN>, UP> ddbetaG;

  // Intermediate variables
  const mat3<T, UP, UP> gammatu;
  const mat3<vec3<T, DN>, UP, UP> dgammatu;
  const vec3<mat3<T, DN, DN>, DN> Gammatl;
  const vec3<vec3<vec3<T, UP>, DN>, DN> Gammatlu;
  const vec3<mat3<T, DN, DN>, UP> Gammat;
  const vec3<T, UP> Gamtd;
  const mat3<T, DN, DN> DDchi;
  const mat3<T, UP, UP> gu;
  const mat3<vec3<T, DN>, DN, DN> dg;
  const vec3<mat3<T, DN, DN>, DN> Gammal;
  const vec3<mat3<T, DN, DN>, UP> Gamma;
  const mat3<T, DN, DN> DDalphaG;
  const mat3<T, DN, DN> Rchi;
  const mat3<T, DN, DN> Rt;
  const mat3<T, DN, DN> R;
  const T Rsc;
  const mat3<T, UP, UP> Atu;
  const mat3<vec3<T, DN>, UP, UP> dAtu;
  const T traceSij;

  // Constraints
  const vec3<T, UP> ZtC;
  const T HC;
  const vec3<T, UP> MtC;
  const T allC;

  // RHS variables
  const T chi_rhs;
  const mat3<T, DN, DN> gammat_rhs;
  const T Kh_rhs;
  const mat3<T, DN, DN> At_rhs;
  const vec3<T, UP> Gamt_rhs;
  const T Theta_rhs;
  const T alphaG_rhs;
  const vec3<T, UP> betaG_rhs;

  // ADM RHS variables
  const mat3<T, DN, DN> K_rhs;
  const T dtalpha_rhs;
  const vec3<T, UP> dtbeta_rhs;

  friend CCTK_ATTRIBUTE_NOINLINE ostream &operator<<(ostream &os,
                                                     const z4c_vars &vars) {
    return os << "z4c_vars{"                               //
              << "kappa1:" << vars.kappa1 << ","           //
              << "kappa2:" << vars.kappa2 << ","           //
              << "f_mu_L:" << vars.f_mu_L << ","           //
              << "f_mu_S:" << vars.f_mu_S << ","           //
              << "eta:" << vars.eta << ","                 //
              << "chi:" << vars.chi << ","                 //
              << "gammat:" << vars.gammat << ","           //
              << "Kh:" << vars.Kh << ","                   //
              << "At:" << vars.At << ","                   //
              << "Gamt:" << vars.Gamt << ","               //
              << "Theta:" << vars.Theta << ","             //
              << "alphaG:" << vars.alphaG << ","           //
              << "betaG:" << vars.betaG << ","             //
              << "eTtt:" << vars.eTtt << ","               //
              << "eTti:" << vars.eTti << ","               //
              << "eTij:" << vars.eTij << ","               //
              << "rho:" << vars.rho << ","                 //
              << "Si:" << vars.Si << ","                   //
              << "Sij:" << vars.Sij << ","                 //
              << "g:" << vars.g << ","                     //
              << "K:" << vars.K << ","                     //
              << "alp:" << vars.alp << ","                 //
              << "dtalp:" << vars.dtalp << ","             //
              << "dtbeta:" << vars.dtbeta << ","           //
              << "dchi:" << vars.dchi << ","               //
              << "ddchi:" << vars.ddchi << ","             //
              << "dgammat:" << vars.dgammat << ","         //
              << "ddgammat:" << vars.ddgammat << ","       //
              << "dKh:" << vars.dKh << ","                 //
              << "dAt:" << vars.dAt << ","                 //
              << "dGamt:" << vars.dGamt << ","             //
              << "dTheta:" << vars.dTheta << ","           //
              << "dalphaG:" << vars.dalphaG << ","         //
              << "ddalphaG:" << vars.ddalphaG << ","       //
              << "dbetaG:" << vars.dbetaG << ","           //
              << "ddbetaG:" << vars.ddbetaG << ","         //
              << "gammatu:" << vars.gammatu << ","         //
              << "dgammatu:" << vars.dgammatu << ","       //
              << "Gammatl:" << vars.Gammatl << ","         //
              << "Gammat:" << vars.Gammat << ","           //
              << "Gamtd:" << vars.Gamtd << ","             //
              << "DDchi:" << vars.DDchi << ","             //
              << "gu:" << vars.gu << ","                   //
              << "dg:" << vars.dg << ","                   //
              << "Gammal:" << vars.Gammal << ","           //
              << "Gamma:" << vars.Gamma << ","             //
              << "DDalphaG:" << vars.DDalphaG << ","       //
              << "Rchi:" << vars.Rchi << ","               //
              << "Rt:" << vars.Rt << ","                   //
              << "R:" << vars.R << ","                     //
              << "Rsc:" << vars.Rsc << ","                 //
              << "Atu:" << vars.Atu << ","                 //
              << "dAtu:" << vars.dAtu << ","               //
              << "traceSij:" << vars.traceSij << ","       //
              << "ZtC:" << vars.ZtC << ","                 //
              << "HC:" << vars.HC << ","                   //
              << "MtC:" << vars.MtC << ","                 //
              << "allC:" << vars.allC << ","               //
              << "chi_rhs:" << vars.chi_rhs << ","         //
              << "gammat_rhs:" << vars.gammat_rhs << ","   //
              << "Kh_rhs:" << vars.Kh_rhs << ","           //
              << "At_rhs:" << vars.At_rhs << ","           //
              << "Gamt_rhs:" << vars.Gamt_rhs << ","       //
              << "Theta_rhs:" << vars.Theta_rhs << ","     //
              << "alphaG_rhs:" << vars.alphaG_rhs << ","   //
              << "betaG_rhs:" << vars.betaG_rhs << ","     //
              << "K_rhs:" << vars.K_rhs << ","             //
              << "dtalpha_rhs:" << vars.dtalpha_rhs << "," //
              << "dtbeta_rhs:" << vars.dtbeta_rhs << ","   //
              << "}";
  }

  // See arXiv:1212.2901 [gr-qc]
  Z4C_INLINE Z4C_GPU z4c_vars(
      const T &kappa1, const T &kappa2, const T &f_mu_L, const T &f_mu_S,
      const T &eta,
      //
      const T &chi, const vec3<T, DN> &dchi, const mat3<T, DN, DN> &ddchi, //
      const mat3<T, DN, DN> &gammat, const mat3<vec3<T, DN>, DN, DN> &dgammat,
      const mat3<mat3<T, DN, DN>, DN, DN> &ddgammat,                   //
      const T &Kh, const vec3<T, DN> &dKh,                             //
      const mat3<T, DN, DN> &At, const mat3<vec3<T, DN>, DN, DN> &dAt, //
      const vec3<T, UP> &Gamt, const vec3<vec3<T, DN>, UP> &dGamt,     //
      const T &Theta, const vec3<T, DN> &dTheta,                       //
      const T &alphaG, const vec3<T, DN> &dalphaG,
      const mat3<T, DN, DN> &ddalphaG, //
      const vec3<T, UP> &betaG, const vec3<vec3<T, DN>, UP> &dbetaG,
      const vec3<mat3<T, DN, DN>, UP> &ddbetaG,
      //
      const T &eTtt, const vec3<T, DN> &eTti, const mat3<T, DN, DN> &eTij)
      : z4c_vars_noderivs<T>(kappa1, kappa2, f_mu_L, f_mu_S, eta, //
                             chi, gammat, Kh, At, Gamt, Theta, alphaG, betaG,
                             eTtt, eTti, eTij),
        // Derivatives of Z4c variables
        dchi(dchi), ddchi(ddchi),             //
        dgammat(dgammat), ddgammat(ddgammat), //
        dKh(dKh),                             //
        dAt(dAt),                             //
        dGamt(dGamt),                         //
        dTheta(dTheta),                       //
        dalphaG(dalphaG), ddalphaG(ddalphaG), //
        dbetaG(dbetaG), ddbetaG(ddbetaG),
        // Intermediate variables
        gammatu(gammat.inv(1)),                   //
        dgammatu(calc_dgu(gammatu, dgammat)),     //
        Gammatl(calc_gammal(dgammat)),            //
        Gammatlu(calc_gammalu(gammatu, Gammatl)), //
        Gammat(calc_gamma(gammatu, Gammatl)),     //
        Gamtd([&] Z4C_INLINE(int a) {
          return sum2([&] Z4C_INLINE(int x, int y) {
            return gammatu(x, y) * Gammat(a)(x, y);
          });
        }), //
        DDchi([&] Z4C_INLINE(int a, int b) {
          return ddchi(a, b) //
                 - sum1([&] Z4C_INLINE(int x) {
                     return Gammat(x)(a, b) * dchi(x);
                   });
        }),                                                               //
        gu([&] Z4C_INLINE(int a, int b) { return chi * gammatu(a, b); }), //
        dg([&] Z4C_INLINE(int a, int b) {
          return vec3<T, DN>([&] Z4C_INLINE(int c) {
            return -dchi(c) / pow2(chi) * gammat(a, b) //
                   + 1 / chi * dgammat(a, b)(c);
          });
        }),                            //
        Gammal(calc_gammal(dg)),       //
        Gamma(calc_gamma(gu, Gammal)), //
        DDalphaG([&] Z4C_INLINE(int a, int b) {
          return ddalphaG(a, b) //
                 - sum1([&] Z4C_INLINE(int x) {
                     return Gamma(x)(a, b) * dalphaG(x);
                   });
          // return ddalphaG(a, b) //
          //        - sum1([&]Z4C_INLINE(int x) { return Gammat(x)(a, b) *
          //        dalphaG(x); })
          //        //
          //        // chipsipower = -4
          //        // oochipsipower = -1
          //        // df = oochipsipower dchi(a) / chi
          //        // ddf == oochipsipower ddchi / chi - (-4) df df
          //        + 1 / (2 * chi) *
          //              (dchi(a) * dalphaG(b)    //
          //               + dchi(b) * dalphaG(a)) //
          //        + 1 / (4 * chi) * sum2([&]Z4C_INLINE(int x, int y) {
          //            return gammat(a, b) * gammatu(x, y) * dchi(x) *
          //            dalphaG(y);
          //          });
        }),
        // (8)
        Rchi([&] Z4C_INLINE(int a, int b) {
          return 1 / (2 * chi) * DDchi(a, b) //
                 + 1 / (2 * chi) * gammat(a, b) *
                       sum2([&] Z4C_INLINE(int x, int y) {
                         return gammatu(x, y) * DDchi(x, y);
                       })                                  //
                 - 1 / (4 * pow2(chi)) * dchi(a) * dchi(b) //
                 - 3 / (4 * pow2(chi)) * gammat(a, b) *
                       sum2([&] Z4C_INLINE(int x, int y) {
                         return gammatu(x, y) * dchi(x) * dchi(y);
                       });
        }),
        // (9)
        Rt([&] Z4C_INLINE(int a, int b) {
          return //
              -1 / T(2) * sum2([&] Z4C_INLINE(int x, int y) {
                return gammatu(x, y) * ddgammat(a, b)(x, y);
              }) //
              + 1 / T(2) * sum1([&] Z4C_INLINE(int x) {
                  return (gammat(x, a) * dGamt(x)(b) //
                          + gammat(x, b) * dGamt(x)(a));
                }) //
              + 1 / T(2) * sum1([&] Z4C_INLINE(int x) {
                  return (Gamtd(x) * Gammatl(a)(b, x) //
                          + Gamtd(x) * Gammatl(b)(a, x));
                }) //
              // + sum2([&] Z4C_INLINE(int y, int z) {
              //     return gammatu(y, z) * sum1([&] Z4C_INLINE(int x) {
              //              return (Gammat(x)(a, y) * Gammatl(b)(x, z)   //
              //                      + Gammat(x)(b, y) * Gammatl(a)(x, z) //
              //                      + Gammat(x)(a, y) * Gammatl(x)(b, z));
              //            });
              //   }) //
              + sum2([&] Z4C_INLINE(int x, int y) {
                  return Gammat(x)(a, y) * Gammatlu(b)(x)(y)   //
                         + Gammat(x)(b, y) * Gammatlu(a)(x)(y) //
                         + Gammat(x)(a, y) * Gammatlu(x)(b)(y);
                }) //
              ;
        }),
        // (7)
        R([&] Z4C_INLINE(int a, int b) { return Rchi(a, b) + Rt(a, b); }), //
        Rsc(R.trace(gu)),                                                  //
        Atu([&] Z4C_INLINE(int a, int b) {
          return sum1([&] Z4C_INLINE(int x) {
            return gammatu(a, x) * sum1([&] Z4C_INLINE(int y) {
                     return gammatu(b, y) * At(x, y);
                   });
          });
        }),                                         //
        dAtu(calc_dAu(gammatu, dgammatu, At, dAt)), //
        traceSij(Sij.trace(gu)),
        // Constraints
        // (13)
        ZtC([&] Z4C_INLINE(int a) { return (Gamt(a) - Gamtd(a)) / 2; }), //
        // (14)
        HC(Rsc //
           + sum2([&] Z4C_INLINE(int x, int y) {
               return At(x, y) * Atu(x, y);
             })                              //
           - 2 / T(3) * pow2(Kh + 2 * Theta) //
           - 16 * M_PI * rho),
        // (15)
        MtC([&] Z4C_INLINE(int a) {
          return sum1([&] Z4C_INLINE(int x) { return dAtu(a, x)(x); }) //
                 + sum2([&] Z4C_INLINE(int x, int y) {
                     return Gammat(a)(x, y) * Atu(x, y);
                   }) //
                 - 2 / T(3) * sum1([&] Z4C_INLINE(int x) {
                     return gammatu(a, x) * (dKh(x) + 2 * dTheta(x));
                   }) //
                 - 2 / T(3) * sum1([&] Z4C_INLINE(int x) {
                     return Atu(a, x) * dchi(x) / chi;
                   }) //
                 - 8 * M_PI * sum1([&] Z4C_INLINE(int x) {
                     return gammatu(a, x) * Si(x);
                   });
        }),
        // arXiv:1111.2177, (73)
        allC(
            sqrt(fmax(T(0), pow2(HC) //
                                + sum1([&] Z4C_INLINE(int x) {
                                    return MtC(x) * sum1([&] Z4C_INLINE(int y) {
                                             return gammat(x, y) * MtC(y);
                                           });
                                  })          //
                                + pow2(Theta) //
                                + 2 * sum1([&] Z4C_INLINE(int x) {
                                    return ZtC(x) * sum1([&] Z4C_INLINE(int y) {
                                             return gammat(x, y) * ZtC(y);
                                           });
                                  })))),
        // RHS
        // (1)
        chi_rhs(2 / T(3) * chi *
                (alphaG * (Kh + 2 * Theta)
                 // chi = detg^(-1/3)
                 // chi^(-3) = detg
                 // chi^(-3/2) = sqrt detg
                 // D_i beta^i = 1/sqrt detg d_i sqrt detg beta^i
                 //            = chi^(3/2) d_i chi^(-3/2) beta^i
                 //            = chi^(3/2) (-3/2 chi^(-5/2) chi,i beta^i
                 //                         + chi^(-3/2) beta^i,i)
                 //            = beta^i,i - 3/2 1/chi beta^i chi,i
                 - sum1([&] Z4C_INLINE(int x) { return dbetaG(x)(x); }))),
        // (2)
        gammat_rhs([&] Z4C_INLINE(int a, int b) {
          return -2 * alphaG * At(a, b) //
                 + sum1([&] Z4C_INLINE(int x) {
                     return gammat(x, a) * dbetaG(x)(b) //
                            + gammat(x, b) * dbetaG(x)(a);
                   }) //
                 - 2 / T(3) * sum1([&] Z4C_INLINE(int x) {
                     return gammat(a, b) * dbetaG(x)(x);
                   });
        }),
        // (3)
        Kh_rhs(sum2([&] Z4C_INLINE(int x, int y) {
                 return -gu(x, y) * DDalphaG(x, y);
               }) //
               + alphaG * (sum2([&] Z4C_INLINE(int x, int y) {
                             return At(x, y) * Atu(x, y);
                           })                                 //
                           + 1 / T(3) * pow2(Kh + 2 * Theta)) //
               + 4 * M_PI * alphaG * (traceSij + rho)         //
               + alphaG * kappa1 * (1 - kappa2) * Theta),
        // (4)
        At_rhs([&] Z4C_INLINE(int a, int b) {
          return chi *
                     ((-DDalphaG(a, b)                      //
                       + alphaG * (R(a, b)                  //
                                   - 8 * M_PI * Sij(a, b))) //
                      - 1 / T(3) * g(a, b) * sum2([&] Z4C_INLINE(int x, int y) {
                          return gu(x, y) *
                                 (-DDalphaG(x, y)     //
                                  + alphaG * (R(x, y) //
                                              - 8 * M_PI * Sij(x, y)));
                        }))                              //
                 + alphaG * ((Kh + 2 * Theta) * At(a, b) //
                             - 2 * sum1([&] Z4C_INLINE(int x) {
                                 return At(x, a) * sum1([&] Z4C_INLINE(int y) {
                                          return gammatu(x, y) * At(y, b);
                                        });
                               })) //
                 + sum1([&] Z4C_INLINE(int x) {
                     return At(x, a) * dbetaG(x)(b) //
                            + At(x, b) * dbetaG(x)(a);
                   }) //
                 - 2 / T(3) * At(a, b) *
                       sum1([&] Z4C_INLINE(int x) { return dbetaG(x)(x); }); //
        }),
        // (5)
        Gamt_rhs([&] Z4C_INLINE(int a) {
          return //
              -2 * sum1([&] Z4C_INLINE(int x) {
                return Atu(a, x) * dalphaG(x);
              }) //
              + 2 * alphaG *
                    (sum2([&] Z4C_INLINE(int x, int y) {
                       return Gammat(a)(x, y) * Atu(x, y);
                     }) //
                     - 3 / T(2) / chi * sum1([&] Z4C_INLINE(int x) {
                         return Atu(a, x) * dchi(x);
                       }) //
                     - 1 / T(3) * sum1([&] Z4C_INLINE(int x) {
                         return gammatu(a, x) * (2 * dKh(x) + dTheta(x));
                       }) //
                     - 8 * M_PI * sum1([&] Z4C_INLINE(int x) {
                         return gammatu(a, x) * Si(x);
                       })) //
              + sum2([&] Z4C_INLINE(int x, int y) {
                  return gammatu(x, y) * ddbetaG(a)(x, y);
                }) //
              + 1 / T(3) * sum1([&] Z4C_INLINE(int x) {
                  return gammatu(a, x) * sum1([&] Z4C_INLINE(int y) {
                           return ddbetaG(y)(x, y);
                         });
                }) //
              -
              sum1([&] Z4C_INLINE(int x) { return Gamtd(x) * dbetaG(a)(x); }) //
              + 2 / T(3) * Gamtd(a) *
                    sum1([&] Z4C_INLINE(int x) { return dbetaG(x)(x); }) //
              - 2 * alphaG * kappa1 * (Gamt(a) - Gamtd(a));
        }),
        // (6)
        Theta_rhs(1 / T(2) * alphaG *
                      (Rsc //
                       - sum2([&] Z4C_INLINE(int x, int y) {
                           return At(x, y) * Atu(x, y);
                         })                               //
                       + 2 / T(3) * pow2(Kh + 2 * Theta)) //
                  - alphaG * (8 * M_PI * rho              //
                              + kappa1 * (2 + kappa2) * Theta)),
        //
        alphaG_rhs(dtalp), betaG_rhs(dtbeta),
        //
        K_rhs([&] Z4C_INLINE(int a, int b) {
          return -1 / pow(chi, 2) * chi_rhs *
                     (At(a, b) + (Kh + 2 * Theta) / 3 * gammat(a, b)) +
                 1 / chi *
                     (At_rhs(a, b) +
                      (Kh_rhs + 2 * Theta_rhs) / 3 * gammat(a, b)) +
                 1 / chi * (At(a, b) + (Kh + 2 * Theta) / 3 * gammat_rhs(a, b));
        }),
        //
        dtalpha_rhs([&] Z4C_INLINE {
          const T mu_L = f_mu_L / alphaG;
          return -2 * alphaG * alphaG_rhs * mu_L * Kh //
                 - pow2(alphaG) * mu_L * Kh_rhs;
        }()),
        //
        dtbeta_rhs([&] Z4C_INLINE(int a) {
          const T mu_S = f_mu_S / pow2(alphaG);
          return 2 * alphaG * alphaG_rhs * mu_S * Gamt(a) //
                 + pow2(alphaG) * mu_S * Gamt_rhs(a)      //
                 - eta * betaG_rhs(a);
        })
  //
  {}

  Z4C_INLINE Z4C_GPU z4c_vars(
      const T &kappa1, const T &kappa2, const T &f_mu_L, const T &f_mu_S,
      const T &eta,
      //
      const GF3D2<const T> &gf_chi_,
      //
      const GF3D2<const T> &gf_gammatxx_, const GF3D2<const T> &gf_gammatxy_,
      const GF3D2<const T> &gf_gammatxz_, const GF3D2<const T> &gf_gammatyy_,
      const GF3D2<const T> &gf_gammatyz_, const GF3D2<const T> &gf_gammatzz_,
      //
      const GF3D2<const T> &gf_Kh_,
      //
      const GF3D2<const T> &gf_Atxx_, const GF3D2<const T> &gf_Atxy_,
      const GF3D2<const T> &gf_Atxz_, const GF3D2<const T> &gf_Atyy_,
      const GF3D2<const T> &gf_Atyz_, const GF3D2<const T> &gf_Atzz_,
      //
      const GF3D2<const T> &gf_Gamtx_, const GF3D2<const T> &gf_Gamty_,
      const GF3D2<const T> &gf_Gamtz_,
      //
      const GF3D2<const T> &gf_Theta_,
      //
      const GF3D2<const T> &gf_alphaG_,
      //
      const GF3D2<const T> &gf_betaGx_, const GF3D2<const T> &gf_betaGy_,
      const GF3D2<const T> &gf_betaGz_,
      //
      const GF3D2<const T> &gf_eTtt_,
      //
      const GF3D2<const T> &gf_eTtx_, const GF3D2<const T> &gf_eTty_,
      const GF3D2<const T> &gf_eTtz_,
      //
      const GF3D2<const T> &gf_eTxx_, const GF3D2<const T> &gf_eTxy_,
      const GF3D2<const T> &gf_eTxz_, const GF3D2<const T> &gf_eTyy_,
      const GF3D2<const T> &gf_eTyz_, const GF3D2<const T> &gf_eTzz_,
      //
      const vect<int, 3> &I, const vec3<T, UP> &dx)
      : z4c_vars(kappa1, kappa2, f_mu_L, f_mu_S, eta,
                 //
                 gf_chi_(I), deriv(gf_chi_, I, dx), deriv2(gf_chi_, I, dx),
                 //
                 mat3<T, DN, DN>(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_,
                                 gf_gammatyy_, gf_gammatyz_, gf_gammatzz_, I),
                 mat3<vec3<T, DN>, DN, DN>{
                     deriv(gf_gammatxx_, I, dx),
                     deriv(gf_gammatxy_, I, dx),
                     deriv(gf_gammatxz_, I, dx),
                     deriv(gf_gammatyy_, I, dx),
                     deriv(gf_gammatyz_, I, dx),
                     deriv(gf_gammatzz_, I, dx),
                 },
                 mat3<mat3<T, DN, DN>, DN, DN>{
                     deriv2(gf_gammatxx_, I, dx),
                     deriv2(gf_gammatxy_, I, dx),
                     deriv2(gf_gammatxz_, I, dx),
                     deriv2(gf_gammatyy_, I, dx),
                     deriv2(gf_gammatyz_, I, dx),
                     deriv2(gf_gammatzz_, I, dx),
                 },
                 //
                 gf_Kh_(I), deriv(gf_Kh_, I, dx),
                 //
                 mat3<T, DN, DN>(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_,
                                 gf_Atyz_, gf_Atzz_, I),
                 mat3<vec3<T, DN>, DN, DN>{
                     deriv(gf_Atxx_, I, dx),
                     deriv(gf_Atxy_, I, dx),
                     deriv(gf_Atxz_, I, dx),
                     deriv(gf_Atyy_, I, dx),
                     deriv(gf_Atyz_, I, dx),
                     deriv(gf_Atzz_, I, dx),
                 },
                 //
                 vec3<T, UP>(gf_Gamtx_, gf_Gamty_, gf_Gamtz_, I),
                 vec3<vec3<T, DN>, UP>{
                     deriv(gf_Gamtx_, I, dx),
                     deriv(gf_Gamty_, I, dx),
                     deriv(gf_Gamtz_, I, dx),
                 },
                 //
                 gf_Theta_(I), deriv(gf_Theta_, I, dx),
                 //
                 gf_alphaG_(I), deriv(gf_alphaG_, I, dx),
                 deriv2(gf_alphaG_, I, dx),
                 //
                 vec3<T, UP>(gf_betaGx_, gf_betaGy_, gf_betaGz_, I),
                 vec3<vec3<T, DN>, UP>{
                     deriv(gf_betaGx_, I, dx),
                     deriv(gf_betaGy_, I, dx),
                     deriv(gf_betaGz_, I, dx),
                 },
                 vec3<mat3<T, DN, DN>, UP>{
                     deriv2(gf_betaGx_, I, dx),
                     deriv2(gf_betaGy_, I, dx),
                     deriv2(gf_betaGz_, I, dx),
                 },
                 //
                 gf_eTtt_(I),
                 //
                 vec3<T, DN>(gf_eTtx_, gf_eTty_, gf_eTtz_, I),
                 //
                 mat3<T, DN, DN>(gf_eTxx_, gf_eTxy_, gf_eTxz_, gf_eTyy_,
                                 gf_eTyz_, gf_eTzz_, I))
  //
  {}

  template <int NI, int NJ, int NK>
  Z4C_INLINE Z4C_GPU z4c_vars(const T &kappa1, const T &kappa2, const T &f_mu_L,
                              const T &f_mu_S, const T &eta,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_chi_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_gammatxx_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_gammatxy_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_gammatxz_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_gammatyy_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_gammatyz_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_gammatzz_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Kh_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Atxx_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Atxy_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Atxz_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Atyy_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Atyz_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Atzz_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Gamtx_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Gamty_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Gamtz_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_Theta_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_alphaG_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_betaGx_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_betaGy_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_betaGz_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTtt_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTtx_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTty_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTtz_,
                              //
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTxx_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTxy_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTxz_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTyy_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTyz_,
                              const GF3D3ptr<const T, NI, NJ, NK> &gf_eTzz_,
                              //
                              const vect<int, 3> &I, const vec3<T, UP> &dx)
      : z4c_vars(kappa1, kappa2, f_mu_L, f_mu_S, eta,
                 //
                 gf_chi_(I), deriv(gf_chi_, I, dx), deriv2(gf_chi_, I, dx),
                 //
                 mat3<T, DN, DN>(gf_gammatxx_, gf_gammatxy_, gf_gammatxz_,
                                 gf_gammatyy_, gf_gammatyz_, gf_gammatzz_, I),
                 mat3<vec3<T, DN>, DN, DN>{
                     deriv(gf_gammatxx_, I, dx),
                     deriv(gf_gammatxy_, I, dx),
                     deriv(gf_gammatxz_, I, dx),
                     deriv(gf_gammatyy_, I, dx),
                     deriv(gf_gammatyz_, I, dx),
                     deriv(gf_gammatzz_, I, dx),
                 },
                 mat3<mat3<T, DN, DN>, DN, DN>{
                     deriv2(gf_gammatxx_, I, dx),
                     deriv2(gf_gammatxy_, I, dx),
                     deriv2(gf_gammatxz_, I, dx),
                     deriv2(gf_gammatyy_, I, dx),
                     deriv2(gf_gammatyz_, I, dx),
                     deriv2(gf_gammatzz_, I, dx),
                 },
                 //
                 gf_Kh_(I), deriv(gf_Kh_, I, dx),
                 //
                 mat3<T, DN, DN>(gf_Atxx_, gf_Atxy_, gf_Atxz_, gf_Atyy_,
                                 gf_Atyz_, gf_Atzz_, I),
                 mat3<vec3<T, DN>, DN, DN>{
                     deriv(gf_Atxx_, I, dx),
                     deriv(gf_Atxy_, I, dx),
                     deriv(gf_Atxz_, I, dx),
                     deriv(gf_Atyy_, I, dx),
                     deriv(gf_Atyz_, I, dx),
                     deriv(gf_Atzz_, I, dx),
                 },
                 //
                 vec3<T, UP>(gf_Gamtx_, gf_Gamty_, gf_Gamtz_, I),
                 vec3<vec3<T, DN>, UP>{
                     deriv(gf_Gamtx_, I, dx),
                     deriv(gf_Gamty_, I, dx),
                     deriv(gf_Gamtz_, I, dx),
                 },
                 //
                 gf_Theta_(I), deriv(gf_Theta_, I, dx),
                 //
                 gf_alphaG_(I), deriv(gf_alphaG_, I, dx),
                 deriv2(gf_alphaG_, I, dx),
                 //
                 vec3<T, UP>(gf_betaGx_, gf_betaGy_, gf_betaGz_, I),
                 vec3<vec3<T, DN>, UP>{
                     deriv(gf_betaGx_, I, dx),
                     deriv(gf_betaGy_, I, dx),
                     deriv(gf_betaGz_, I, dx),
                 },
                 vec3<mat3<T, DN, DN>, UP>{
                     deriv2(gf_betaGx_, I, dx),
                     deriv2(gf_betaGy_, I, dx),
                     deriv2(gf_betaGz_, I, dx),
                 },
                 //
                 gf_eTtt_(I),
                 //
                 vec3<T, DN>(gf_eTtx_, gf_eTty_, gf_eTtz_, I),
                 //
                 mat3<T, DN, DN>(gf_eTxx_, gf_eTxy_, gf_eTxz_, gf_eTyy_,
                                 gf_eTyz_, gf_eTzz_, I))
  //
  {}
};

} // namespace Z4c

#endif // #ifndef Z4C_HXX
