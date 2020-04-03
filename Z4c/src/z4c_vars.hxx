#ifndef Z4C_HXX
#define Z4C_HXX

#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"

#include <cmath>

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

  // Hydro variables
  const T rho;
  const vec3<T, DN> Si;
  const mat3<T, DN, DN> Sij;

  // ADM variables
  const mat3<T, DN, DN> g; // W = 0
  const mat3<T, DN, DN> k; // W = 0
  const T alp;             // W = 0
  const vec3<T, UP> beta;  // W = 0

  CCTK_ATTRIBUTE_ALWAYS_INLINE
  z4c_vars_noderivs(const T &kappa1, const T &kappa2, const T &f_mu_L,
                    const T &f_mu_S, const T &eta,
                    //
                    const T &chi, const mat3<T, DN, DN> &gammat, const T &Kh,
                    const mat3<T, DN, DN> &At, const vec3<T, UP> &Gamt,
                    const T &Theta, const T &alphaG, const vec3<T, UP> &betaG,
                    //
                    const T &rho, const vec3<T, DN> &Si,
                    const mat3<T, DN, DN> &Sij)
      : kappa1(kappa1), kappa2(kappa2), f_mu_L(f_mu_L), f_mu_S(f_mu_S),
        eta(eta),
        //
        chi(chi), gammat(gammat), Kh(Kh), At(At), Gamt(Gamt), Theta(Theta),
        alphaG(alphaG), betaG(betaG),
        //
        rho(rho), Si(Si), Sij(Sij),
        // ADM variables
        g([&](int a, int b) { return 1 / chi * gammat(a, b); }), //
        k([&](int a, int b) {
          return 1 / chi * (At(a, b) + (Kh + 2 * Theta) / 3 * gammat(a, b));
        }),          //
        alp(alphaG), //
        beta([&](int a) { return betaG(a); })
  //
  {}

  CCTK_ATTRIBUTE_ALWAYS_INLINE
  z4c_vars_noderivs(const T &kappa1, const T &kappa2, const T &f_mu_L,
                    const T &f_mu_SL, const T &eta,
                    //
                    const GF3D<const T, 0, 0, 0> &gf_chi_,
                    //
                    const GF3D<const T, 0, 0, 0> &gf_gammatxx_,
                    const GF3D<const T, 0, 0, 0> &gf_gammatxy_,
                    const GF3D<const T, 0, 0, 0> &gf_gammatxz_,
                    const GF3D<const T, 0, 0, 0> &gf_gammatyy_,
                    const GF3D<const T, 0, 0, 0> &gf_gammatyz_,
                    const GF3D<const T, 0, 0, 0> &gf_gammatzz_,
                    //
                    const GF3D<const T, 0, 0, 0> &gf_Kh_,
                    //
                    const GF3D<const T, 0, 0, 0> &gf_Atxx_,
                    const GF3D<const T, 0, 0, 0> &gf_Atxy_,
                    const GF3D<const T, 0, 0, 0> &gf_Atxz_,
                    const GF3D<const T, 0, 0, 0> &gf_Atyy_,
                    const GF3D<const T, 0, 0, 0> &gf_Atyz_,
                    const GF3D<const T, 0, 0, 0> &gf_Atzz_,
                    //
                    const GF3D<const T, 0, 0, 0> &gf_Gamtx_,
                    const GF3D<const T, 0, 0, 0> &gf_Gamty_,
                    const GF3D<const T, 0, 0, 0> &gf_Gamtz_,
                    //
                    const GF3D<const T, 0, 0, 0> &gf_Theta_,
                    //
                    const GF3D<const T, 0, 0, 0> &gf_alphaG_,
                    //
                    const GF3D<const T, 0, 0, 0> &gf_betaGx_,
                    const GF3D<const T, 0, 0, 0> &gf_betaGy_,
                    const GF3D<const T, 0, 0, 0> &gf_betaGz_,
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
                             0,
                             //
                             vec3<T, DN>{0, 0, 0},
                             //
                             mat3<T, DN, DN>{0, 0, 0, 0, 0, 0})
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

  // Hydro variables
  using z4c_vars_noderivs<T>::rho;
  using z4c_vars_noderivs<T>::Si;
  using z4c_vars_noderivs<T>::Sij;

  // ADM variables
  using z4c_vars_noderivs<T>::g;
  using z4c_vars_noderivs<T>::k;
  using z4c_vars_noderivs<T>::alp;
  using z4c_vars_noderivs<T>::beta;

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

  // See arXiv:1212.2901 [gr-qc]
  CCTK_ATTRIBUTE_ALWAYS_INLINE
  z4c_vars(const T &kappa1, const T &kappa2, const T &f_mu_L, const T &f_mu_S,
           const T &eta,
           //
           const T &chi, const vec3<T, DN> &dchi,
           const mat3<T, DN, DN> &ddchi, //
           const mat3<T, DN, DN> &gammat,
           const mat3<vec3<T, DN>, DN, DN> &dgammat,
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
           const T &rho, const vec3<T, DN> &Si, const mat3<T, DN, DN> &Sij)
      : z4c_vars_noderivs<T>(kappa1, kappa2, f_mu_L, f_mu_S, eta, //
                             chi, gammat, Kh, At, Gamt, Theta, alphaG, betaG,
                             rho, Si, Sij),
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
        gammatu(gammat.inv(1)),               //
        dgammatu(calc_dgu(gammatu, dgammat)), //
        Gammatl(calc_gammal(dgammat)),        //
        Gammat(calc_gamma(gammatu, Gammatl)), //
        Gamtd([&](int a) {
          return sum2(
              [&](int x, int y) { return gammatu(x, y) * Gammat(a)(x, y); });
        }), //
        DDchi([&](int a, int b) {
          return ddchi(a, b) //
                 - sum1([&](int x) { return Gammat(x)(a, b) * dchi(x); });
        }),                                                    //
        gu([&](int a, int b) { return chi * gammatu(a, b); }), //
        dg([&](int a, int b) {
          return vec3<T, DN>([&](int c) {
            return -dchi(c) / pow2(chi) * gammat(a, b) //
                   + 1 / chi * dgammat(a, b)(c);
          });
        }),                            //
        Gammal(calc_gammal(dg)),       //
        Gamma(calc_gamma(gu, Gammal)), //
        DDalphaG([&](int a, int b) {
          return ddalphaG(a, b) //
                 - sum1([&](int x) { return Gamma(x)(a, b) * dalphaG(x); });
          // return ddalphaG(a, b) //
          //        - sum1([&](int x) { return Gammat(x)(a, b) * dalphaG(x); })
          //        //
          //        // chipsipower = -4
          //        // oochipsipower = -1
          //        // df = oochipsipower dchi(a) / chi
          //        // ddf == oochipsipower ddchi / chi - (-4) df df
          //        + 1 / (2 * chi) *
          //              (dchi(a) * dalphaG(b)    //
          //               + dchi(b) * dalphaG(a)) //
          //        + 1 / (4 * chi) * sum2([&](int x, int y) {
          //            return gammat(a, b) * gammatu(x, y) * dchi(x) *
          //            dalphaG(y);
          //          });
        }),
        // (8)
        Rchi([&](int a, int b) {
          return 1 / (2 * chi) * DDchi(a, b) //
                 + 1 / (2 * chi) * sum2([&](int x, int y) {
                     return gammat(a, b) * gammatu(x, y) * DDchi(x, y);
                   })                                      //
                 - 1 / (4 * pow2(chi)) * dchi(a) * dchi(b) //
                 - 3 / (4 * pow2(chi)) * sum2([&](int x, int y) {
                     return gammat(a, b) * gammatu(x, y) * dchi(x) * dchi(y);
                   });
        }),
        // (9)
        Rt([&](int a, int b) {
          return sum2([&](int x, int y) {
                   return -1 / T(2) * gammatu(x, y) * ddgammat(a, b)(x, y);
                 }) //
                 + sum1([&](int x) {
                     return 1 / T(2) *
                            (gammat(x, a) * dGamt(x)(b) //
                             + gammat(x, b) * dGamt(x)(a));
                   }) //
                 + sum1([&](int x) {
                     return 1 / T(2) *
                            (Gamtd(x) * Gammatl(a)(b, x) //
                             + Gamtd(x) * Gammatl(b)(a, x));
                   }) //
                 + sum3([&](int x, int y, int z) {
                     return gammatu(y, z) *
                            (Gammat(x)(a, y) * Gammatl(b)(x, z)   //
                             + Gammat(x)(b, y) * Gammatl(a)(x, z) //
                             + Gammat(x)(a, y) * Gammatl(x)(b, z));
                   });
        }),
        // (7)
        R([&](int a, int b) { return Rchi(a, b) + Rt(a, b); }), //
        Rsc(R.trace(gu)),                                       //
        Atu([&](int a, int b) {
          return sum2([&](int x, int y) {
            return gammatu(a, x) * gammatu(b, y) * At(x, y);
          });
        }),                                         //
        dAtu(calc_dAu(gammatu, dgammatu, At, dAt)), //
        traceSij(Sij.trace(gu)),
        // Constraints
        // (13)
        ZtC([&](int a) { return (Gamt(a) - Gamtd(a)) / 2; }), //
        // (14)
        HC(Rsc                                                        //
           + sum2([&](int x, int y) { return At(x, y) * Atu(x, y); }) //
           - 2 / T(3) * pow2(Kh + 2 * Theta)                          //
           - 16 * M_PI * rho),
        // (15)
        MtC([&](int a) {
          return sum1([&](int x) { return dAtu(a, x)(x); }) //
                 + sum2([&](int x, int y) {
                     return Gammat(a)(x, y) * Atu(x, y);
                   }) //
                 - sum1([&](int x) {
                     return 2 / T(3) * gammatu(a, x) * (dKh(x) + 2 * dTheta(x));
                   }) //
                 - sum1([&](int x) {
                     return 2 / T(3) * Atu(a, x) * dchi(x) / chi;
                   }) //
                 -
                 sum1([&](int x) { return 8 * M_PI * gammatu(a, x) * Si(x); });
        }),
        // arXiv:1111.2177, (73)
        allC(sqrt(pow2(HC) //
                  + sum2([&](int x, int y) {
                      return gammat(x, y) * MtC(x) * MtC(y);
                    })          //
                  + pow2(Theta) //
                  + sum2([&](int x, int y) {
                      return 2 * gammat(x, y) * ZtC(x) * ZtC(y);
                    }))),
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
                 - sum1([&](int x) { return dbetaG(x)(x); }))),
        // (2)
        gammat_rhs([&](int a, int b) {
          return -2 * alphaG * At(a, b) //
                 + sum1([&](int x) {
                     return gammat(x, a) * dbetaG(x)(b) //
                            + gammat(x, b) * dbetaG(x)(a);
                   }) //
                 - sum1([&](int x) {
                     return 2 / T(3) * gammat(a, b) * dbetaG(x)(x);
                   });
        }),
        // (3)
        Kh_rhs(sum2([&](int x, int y) { return -gu(x, y) * DDalphaG(x, y); }) //
               + alphaG * (sum2([&](int x, int y) {
                             return At(x, y) * Atu(x, y);
                           })                                 //
                           + 1 / T(3) * pow2(Kh + 2 * Theta)) //
               + 4 * M_PI * alphaG * (traceSij + rho)         //
               + alphaG * kappa1 * (1 - kappa2) * Theta),
        // (4)
        At_rhs([&](int a, int b) {
          return chi * ((-DDalphaG(a, b)                      //
                         + alphaG * (R(a, b)                  //
                                     - 8 * M_PI * Sij(a, b))) //
                        - sum2([&](int x, int y) {
                            return 1 / T(3) * g(a, b) * gu(x, y) *
                                   (-DDalphaG(x, y)     //
                                    + alphaG * (R(x, y) //
                                                - 8 * M_PI * Sij(x, y)));
                          }))                            //
                 + alphaG * ((Kh + 2 * Theta) * At(a, b) //
                             - sum2([&](int x, int y) {
                                 return 2 * gammatu(x, y) * At(x, a) * At(y, b);
                               })) //
                 + sum1([&](int x) {
                     return At(x, a) * dbetaG(x)(b) //
                            + At(x, b) * dbetaG(x)(a);
                   }) //
                 - sum1([&](int x) {
                     return 2 / T(3) * At(a, b) * dbetaG(x)(x);
                   }); //
        }),
        // (5)
        Gamt_rhs([&](int a) {
          return sum1([&](int x) { return -2 * Atu(a, x) * dalphaG(x); }) //
                 + 2 * alphaG *
                       (sum2([&](int x, int y) {
                          return Gammat(a)(x, y) * Atu(x, y);
                        }) //
                        - sum1([&](int x) {
                            return 3 / T(2) * Atu(a, x) * dchi(x) / chi;
                          }) //
                        - sum1([&](int x) {
                            return 1 / T(3) * gammatu(a, x) *
                                   (2 * dKh(x) + dTheta(x));
                          }) //
                        - sum1([&](int x) {
                            return 8 * M_PI * gammatu(a, x) * Si(x);
                          })) //
                 + sum2([&](int x, int y) {
                     return gammatu(x, y) * ddbetaG(a)(x, y);
                   }) //
                 + sum2([&](int x, int y) {
                     return 1 / T(3) * gammatu(a, x) * ddbetaG(y)(x, y);
                   })                                                   //
                 - sum1([&](int x) { return Gamtd(x) * dbetaG(a)(x); }) //
                 + sum1([&](int x) {
                     return 2 / T(3) * Gamtd(a) * dbetaG(x)(x);
                   }) //
                 - 2 * alphaG * kappa1 * (Gamt(a) - Gamtd(a));
        }),
        // (6)
        Theta_rhs(
            1 / T(2) * alphaG *
                (Rsc                                                        //
                 - sum2([&](int x, int y) { return At(x, y) * Atu(x, y); }) //
                 + 2 / T(3) * pow2(Kh + 2 * Theta))                         //
            - alphaG * (8 * M_PI * rho                                      //
                        + kappa1 * (2 + kappa2) * Theta)),
        // (11)
        alphaG_rhs([&] {
          const T mu_L = f_mu_L / alphaG;
          return -pow2(alphaG) * mu_L * Kh;
        }()),
        // (12)
        betaG_rhs([&](int a) {
          const T mu_S = f_mu_S / pow2(alphaG);
          return pow2(alphaG) * mu_S * Gamt(a) //
                 - eta * betaG(a);
        })
  //
  {}

  CCTK_ATTRIBUTE_ALWAYS_INLINE
  z4c_vars(const T &kappa1, const T &kappa2, const T &f_mu_L, const T &f_mu_SL,
           const T &eta,
           //
           const GF3D<const T, 0, 0, 0> &gf_chi_,
           //
           const GF3D<const T, 0, 0, 0> &gf_gammatxx_,
           const GF3D<const T, 0, 0, 0> &gf_gammatxy_,
           const GF3D<const T, 0, 0, 0> &gf_gammatxz_,
           const GF3D<const T, 0, 0, 0> &gf_gammatyy_,
           const GF3D<const T, 0, 0, 0> &gf_gammatyz_,
           const GF3D<const T, 0, 0, 0> &gf_gammatzz_,
           //
           const GF3D<const T, 0, 0, 0> &gf_Kh_,
           //
           const GF3D<const T, 0, 0, 0> &gf_Atxx_,
           const GF3D<const T, 0, 0, 0> &gf_Atxy_,
           const GF3D<const T, 0, 0, 0> &gf_Atxz_,
           const GF3D<const T, 0, 0, 0> &gf_Atyy_,
           const GF3D<const T, 0, 0, 0> &gf_Atyz_,
           const GF3D<const T, 0, 0, 0> &gf_Atzz_,
           //
           const GF3D<const T, 0, 0, 0> &gf_Gamtx_,
           const GF3D<const T, 0, 0, 0> &gf_Gamty_,
           const GF3D<const T, 0, 0, 0> &gf_Gamtz_,
           //
           const GF3D<const T, 0, 0, 0> &gf_Theta_,
           //
           const GF3D<const T, 0, 0, 0> &gf_alphaG_,
           //
           const GF3D<const T, 0, 0, 0> &gf_betaGx_,
           const GF3D<const T, 0, 0, 0> &gf_betaGy_,
           const GF3D<const T, 0, 0, 0> &gf_betaGz_,
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
                 0,
                 //
                 vec3<T, DN>{0, 0, 0},
                 //
                 mat3<T, DN, DN>{0, 0, 0, 0, 0, 0})
  //
  {}
};

} // namespace Z4c

#endif // #ifndef Z4C_HXX
