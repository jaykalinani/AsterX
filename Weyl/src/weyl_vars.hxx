#ifndef WEYL_HXX
#define WEYL_HXX

#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"

#include <cmath>

namespace Weyl {

template <typename T> struct weyl_vars_noderivs {

  // ADM variables
  const mat3<T, DN, DN> gamma;
  const T alpha;
  const vec3<T, UP> beta;

  // Intermediate quantities
  const vec3<T, DN> betal;

  // 4-metric
  const mat4<T, DN, DN> g;

  // Inverse 4-metric
  const T detg;
  const mat4<T, UP, UP> gu;

  weyl_vars_noderivs(const mat3<T, DN, DN> &gamma, const T &alpha,
                     const vec3<T, UP> &beta)
      : gamma(gamma), alpha(alpha), beta(beta),
        //
        betal([&](int a) {
          return sum1([&](int x) { return gamma(a, x) * beta(x); });
        }),
        //
        g([&](int a, int b) {
          if (a == 0 && b == 0)
            return alpha;
          if (a == 0)
            return betal(b);
          if (b == 0)
            return betal(a);
          return gamma(a, b);
        }),
        detg(g.det()), //
        gu(g.inv(detg))
  //
  {}

#if 0
  
  weyl_vars_noderivs(const GF3D<const T, 0, 0, 0> &gf_gammaxx_,
                     const GF3D<const T, 0, 0, 0> &gf_gammaxy_,
                     const GF3D<const T, 0, 0, 0> &gf_gammaxz_,
                     const GF3D<const T, 0, 0, 0> &gf_gammayy_,
                     const GF3D<const T, 0, 0, 0> &gf_gammayz_,
                     const GF3D<const T, 0, 0, 0> &gf_gammazz_,
                     //
                     const GF3D<const T, 0, 0, 0> &gf_alpha_,
                     //
                     const GF3D<const T, 0, 0, 0> &gf_betax_,
                     const GF3D<const T, 0, 0, 0> &gf_betay_,
                     const GF3D<const T, 0, 0, 0> &gf_betaz_,
                     //
                     const vect<int, 3> &I)
      : weyl_vars_noderivs<T>(mat3<T, DN, DN>(gf_gammaxx_, gf_gammaxy_,
                                              gf_gammaxz_, gf_gammayy_,
                                              gf_gammayz_, gf_gammazz_, I),
                              //
                              gf_alpha_(I),
                              //
                              vec3<T, UP>(gf_betax_, gf_betay_, gf_betaz_, I))
  // {}
#endif
};

template <typename T> struct weyl_vars : weyl_vars_noderivs<T> {

  // C++ is tedious:

  // ADM variables
  using weyl_vars_noderivs<T>::alpha;
  using weyl_vars_noderivs<T>::beta;
  using weyl_vars_noderivs<T>::gamma;

  // Intermediate quantities
  using weyl_vars_noderivs<T>::betal;

  // 4-metric
  using weyl_vars_noderivs<T>::g;

  // Inverse 4-metric
  using weyl_vars_noderivs<T>::detg;
  using weyl_vars_noderivs<T>::gu;

  // Time derivatives of ADM variables
  const mat3<T, DN, DN> k;
  const T dtalpha;
  const vec3<T, UP> dtbeta;

  // Spatial derivatives of ADM variables
  const mat3<vec3<T, DN>, DN, DN> dgamma;
  const mat3<mat3<T, DN, DN>, DN, DN> ddgamma;
  const vec3<T, DN> dalpha;

  // Second time derivatives of ADM variables
  const mat3<T, DN, DN> dtk;
  const T dt2alpha;
  const vec3<T, UP> dt2beta;

  // Space-time derivatives of ADM variables
  const mat3<vec3<T, DN>, DN, DN> dk;
  const vec3<T, DN> ddtalpha;
  const vec3<vec3<T, DN>, UP> ddtbeta;

  // Second spatial derivatives of ADM variables
  const mat3<T, DN, DN> ddalpha;
  const vec3<vec3<T, DN>, UP> dbeta;
  const vec3<mat3<T, DN, DN>, UP> ddbeta;

  // Intermediate quantities
  const mat3<T, DN, DN> dtgamma;
  const vec3<T, DN> dtbetal;
  const vec3<vec3<T, DN>, DN> dbetal;

  const mat3<T, DN, DN> dt2gamma;
  const mat3<vec3<T, DN>, DN, DN> ddtgamma;
  const vec3<T, DN> dt2betal;
  const vec3<vec3<T, DN>, DN> ddtbetal;
  const vec3<mat3<T, DN, DN>, DN> ddbetal;

  // Derivatives of 4-metric
  const mat4<vec4<T, DN>, DN, DN> dg;
  const mat4<mat4<T, DN, DN>, DN, DN> ddg;
  const mat4<vec4<T, DN>, UP, UP> dgu;

  // Christoffel symbol
  const vec4<mat4<T, DN, DN>, DN> Gammal;
  const vec4<mat4<T, DN, DN>, UP> Gamma;

  const vec4<mat4<vec4<T, DN>, DN, DN>, DN> dGammal;
  // const vec4<mat4<vec4<T, DN>, DN, DN>, UP> dGamma;

  // Riemann, Ricci, Weyl
  // TODO: Use Rm(a,b,c,d) == Rm(c,d,a,b)
  const amat4<amat4<T, DN, DN>, DN, DN> Rm;
  const mat4<T, DN, DN> R;
  const T Rsc;
  const amat4<amat4<T, DN, DN>, DN, DN> C;

  weyl_vars(
      const mat3<T, DN, DN> &gamma, const T &alpha, const vec3<T, UP> &beta,
      //
      const mat3<T, DN, DN> &k, const T &dtalpha, const vec3<T, UP> &dtbeta,
      //
      const mat3<vec3<T, DN>, DN, DN> &dgamma, const vec3<T, DN> &dalpha,
      const vec3<vec3<T, DN>, UP> &dbeta,
      //
      const mat3<T, DN, DN> &dtk, const T &dt2alpha, const vec3<T, UP> &dt2beta,
      //
      const mat3<vec3<T, DN>, DN, DN> &dk, const vec3<T, DN> &ddtalpha,
      const vec3<vec3<T, DN>, UP> &ddtbeta,
      //
      const mat3<mat3<T, DN, DN>, DN, DN> &ddgamma,
      const mat3<T, DN, DN> &ddalpha, const vec3<mat3<T, DN, DN>, UP> &ddbeta)
      : weyl_vars_noderivs<T>(gamma, alpha, beta),
        // Time derivatives of ADM variables
        k(k), dtalpha(dtalpha), dtbeta(dtbeta),
        // Spatial derivatives of ADM variables
        dgamma(dgamma), ddgamma(ddgamma), dalpha(dalpha),
        // Second time derivatives of ADM variables
        dtk(dtk), dt2alpha(dt2alpha), dt2beta(dt2beta),
        // Space-time derivatives of ADM variables
        dk(dk), ddtalpha(ddtalpha), ddtbeta(ddtbeta),
        // Second spatial derivatives of ADM variables
        ddalpha(ddalpha), dbeta(dbeta), ddbeta(ddbeta),
        //
        // dt gamma_ij = -2 alpha K_ij
        //               + gamma_kj beta^k,i + gamma_ik beta^k,j
        //               + beta^k gamma_ij,k
        dtgamma([&](int a, int b) {
          return -2 * alpha * k(a, b)                                     //
                 + sum1([&](int x) { return gamma(x, b) * dbeta(x)(a); }) //
                 + sum1([&](int x) { return gamma(a, x) * dbeta(x)(b); }) //
                 + sum1([&](int x) { return beta(x) * dgamma(a, b)(x); });
        }),
        dtbetal([&](int a) {
          return sum1([&](int x) {
            return dtgamma(a, x) * beta(x) + gamma(a, x) * dtbeta(x);
          });
        }),
        dbetal([&](int a) {
          return vec3<T, DN>([&](int b) {
            return sum1([&](int x) {
              return dgamma(a, x)(b) * beta(x) + gamma(a, x) * dbeta(x)(b);
            });
          });
        }),
        //
        dt2gamma([&](int a, int b) {
          return -2 * dtalpha * k(a, b)  //
                 - 2 * alpha * dtk(a, b) //
                 + sum1([&](int x) {
                     return dtgamma(x, b) * dbeta(x)(a) +
                            gamma(x, b) * ddtbeta(x)(a);
                   }) //
                 + sum1([&](int x) {
                     return dtgamma(a, x) * dbeta(x)(b) +
                            gamma(a, x) * ddtbeta(x)(b);
                   }) //
                 + sum1([&](int x) {
                     return dtbeta(x) * dgamma(a, b)(x) +
                            beta(x) * ddtgamma(a, b)(x);
                   });
        }),
        ddtgamma([&](int a, int b) {
          return vec3<T, DN>([&](int c) {
            return -2 * dalpha(c) * k(a, b)  //
                   - 2 * alpha * dk(a, b)(c) //
                   + sum1([&](int x) {
                       return dgamma(x, b)(c) * dbeta(x)(a) +
                              gamma(x, b) * ddbeta(x)(a, c);
                     }) //
                   + sum1([&](int x) {
                       return dgamma(a, x)(c) * dbeta(x)(b) +
                              gamma(a, x) * ddbeta(x)(b, c);
                     }) //
                   + sum1([&](int x) {
                       return dbeta(x)(c) * dgamma(a, b)(x) +
                              beta(x) * ddgamma(a, b)(x, c);
                     });
          });
        }),
        dt2betal([&](int a) {
          return sum1([&](int x) {
            return dt2gamma(a, x) * beta(x)        //
                   + 2 * dtgamma(a, x) * dtbeta(x) //
                   + dtgamma(a, x) * dt2beta(x);
          });
        }),
        ddtbetal([&](int a) {
          return vec3<T, DN>([&](int b) {
            return sum1([&](int x) {
              return ddtgamma(a, x)(b) * beta(x)   //
                     + dgamma(a, x)(b) * dtbeta(x) //
                     + dtgamma(a, x) * dbeta(x)(b) //
                     + gamma(a, x) * ddtbeta(x)(b);
            });
          });
        }),
        ddbetal([&](int a) {
          return mat3<T, DN, DN>([&](int b, int c) {
            return sum1([&](int x) {
              return ddgamma(a, x)(b, c) * beta(x)   //
                     + dgamma(a, x)(b) * dbeta(x)(c) //
                     + dgamma(a, x)(c) * dbeta(x)(b) //
                     + gamma(a, x) * ddbeta(x)(b, c);
            });
          });
        }),
        //
        dg([&](int a, int b) {
          return vec4<T, DN>([&](int c) {
            if (c == 0) {
              if (a == 0 && b == 0)
                return dtalpha;
              if (a == 0)
                return dtbetal(b - 1);
              if (b == 0)
                return dtbetal(a - 1);
              return dtgamma(a - 1, b - 1);
            }
            if (a == 0 && b == 0)
              return dalpha(c - 1);
            if (a == 0)
              return dbetal(b - 1)(c - 1);
            if (b == 0)
              return dbetal(a - 1)(c - 1);
            return dgamma(a - 1, b - 1)(c - 1);
          });
        }),
        //
        ddg([&](int a, int b) {
          return mat4<T, DN, DN>([&](int c, int d) {
            if (c == 0 && d == 0) {
              if (a == 0 && b == 0)
                return dt2alpha;
              if (a == 0)
                return dt2betal(b - 1);
              if (b == 0)
                return dt2betal(a - 1);
              return dt2gamma(a - 1, b - 1);
            }
            if (c == 0) {
              if (a == 0 && b == 0)
                return ddtalpha(d - 1);
              if (a == 0)
                return ddtbetal(b - 1)(d - 1);
              if (b == 0)
                return ddtbetal(a - 1)(d - 1);
              return ddtgamma(a - 1, b - 1)(d - 1);
            }
            if (d == 0) {
              if (a == 0 && b == 0)
                return ddtalpha(c - 1);
              if (a == 0)
                return ddtbetal(b - 1)(c - 1);
              if (b == 0)
                return ddtbetal(a - 1)(c - 1);
              return ddtgamma(a - 1, b - 1)(c - 1);
            }
            if (a == 0 && b == 0)
              return ddalpha(c - 1, d - 1);
            if (a == 0)
              return ddbetal(b - 1)(c - 1, d - 1);
            if (b == 0)
              return ddbetal(a - 1)(c - 1, d - 1);
            return ddgamma(a - 1, b - 1)(c - 1, d - 1);
          });
        }),
        //
        dgu(calc_dgu(gu, dg)),
        //
        Gammal(calc_gammal(dg)),       //
        Gamma(calc_gamma(gu, Gammal)), //
        //
        dGammal(calc_dgammal(ddg)), //
        // dGamma(calc_dgamma(gu,dgu,Gammal,dGammal)),   //
        //
        Rm(calc_riemann(Gammal, Gamma, dGammal)), //
        R(calc_ricci(gu, Rm)),                    //
        Rsc(R.trace(gu)),                         //
        C(calc_weyl(g, Rm, R, Rsc))
  //
  {}

#if 0
  weyl_vars(const GF3D<const T, 0, 0, 0> &gf_gammaxx_,
            const GF3D<const T, 0, 0, 0> &gf_gammaxy_,
            const GF3D<const T, 0, 0, 0> &gf_gammaxz_,
            const GF3D<const T, 0, 0, 0> &gf_gammayy_,
            const GF3D<const T, 0, 0, 0> &gf_gammayz_,
            const GF3D<const T, 0, 0, 0> &gf_gammazz_,
            //
            const GF3D<const T, 0, 0, 0> &gf_alpha_,
            //
            const GF3D<const T, 0, 0, 0> &gf_betax_,
            const GF3D<const T, 0, 0, 0> &gf_betay_,
            const GF3D<const T, 0, 0, 0> &gf_betaz_,
            //
            const GF3D<const T, 0, 0, 0> &gf_kxx_,
            const GF3D<const T, 0, 0, 0> &gf_kxy_,
            const GF3D<const T, 0, 0, 0> &gf_kxz_,
            const GF3D<const T, 0, 0, 0> &gf_kyy_,
            const GF3D<const T, 0, 0, 0> &gf_kyz_,
            const GF3D<const T, 0, 0, 0> &gf_kzz_,
            //
            const GF3D<const T, 0, 0, 0> &gf_dtalpha_,
            //
            const GF3D<const T, 0, 0, 0> &gf_dtbetax_,
            const GF3D<const T, 0, 0, 0> &gf_dtbetay_,
            const GF3D<const T, 0, 0, 0> &gf_dtbetaz_,
            //
            const GF3D<const T, 0, 0, 0> &gf_dtkxx_,
            const GF3D<const T, 0, 0, 0> &gf_dtkxy_,
            const GF3D<const T, 0, 0, 0> &gf_dtkxz_,
            const GF3D<const T, 0, 0, 0> &gf_dtkyy_,
            const GF3D<const T, 0, 0, 0> &gf_dtkyz_,
            const GF3D<const T, 0, 0, 0> &gf_dtkzz_,
            //
            const GF3D<const T, 0, 0, 0> &gf_dt2alpha_,
            //
            const GF3D<const T, 0, 0, 0> &gf_dt2betax_,
            const GF3D<const T, 0, 0, 0> &gf_dt2betay_,
            const GF3D<const T, 0, 0, 0> &gf_dt2betaz_,
            //
            const vect<int, 3> &I, const vec3<T, UP> &dx)
      : weyl_vars(mat3<T, DN, DN>(gf_gammaxx_, gf_gammaxy_, gf_gammaxz_,
                                  gf_gammayy_, gf_gammayz_, gf_gammazz_, I),
                  //
                  gf_alpha_(I),
                  //
                  vec3<T, UP>(gf_betax_, gf_betay_, gf_betaz_, I),
                  //
                  mat3<T, DN, DN>(gf_kxx_, gf_kxy_, gf_kxz_, gf_kyy_, gf_kyz_,
                                  gf_kzz_, I),
                  //
                  gf_dtalpha_(I),
                  //
                  vec3<T, UP>(gf_dtbetax_, gf_dtbetay_, gf_dtbetaz_, I),
                  //
                  mat3<vec3<T, DN>, DN, DN>{
                      deriv(gf_gammaxx_, I, dx),
                      deriv(gf_gammaxy_, I, dx),
                      deriv(gf_gammaxz_, I, dx),
                      deriv(gf_gammayy_, I, dx),
                      deriv(gf_gammayz_, I, dx),
                      deriv(gf_gammazz_, I, dx),
                  },
                  //
                  deriv(gf_alpha_, I, dx),
                  //
                  vec3<vec3<T, DN>, UP>{
                      deriv(gf_betax_, I, dx),
                      deriv(gf_betay_, I, dx),
                      deriv(gf_betaz_, I, dx),
                  },
                  //
                  mat3<T, DN, DN>(gf_dtkxx_, gf_dtkxy_, gf_dtkxz_, gf_dtkyy_,
                                  gf_dtkyz_, gf_dtkzz_, I),
                  //
                  gf_dt2alpha_(I),
                  //
                  vec3<T, UP>(gf_dt2betax_, gf_dt2betay_, gf_dt2betaz_, I),
                  //
                  mat3<vec3<T, DN>, DN, DN>{
                      deriv(gf_kxx_, I, dx),
                      deriv(gf_kxy_, I, dx),
                      deriv(gf_kxz_, I, dx),
                      deriv(gf_kyy_, I, dx),
                      deriv(gf_kyz_, I, dx),
                      deriv(gf_kzz_, I, dx),
                  },
                  //
                  deriv(gf_dtalpha_, I, dx),
                  //
                  vec3<vec3<T, DN>, UP>{
                      deriv(gf_dtbetax_, I, dx),
                      deriv(gf_dtbetay_, I, dx),
                      deriv(gf_dtbetaz_, I, dx),
                  },
                  //
                  mat3<mat3<T, DN, DN>, DN, DN>{
                      deriv2(gf_gammaxx_, I, dx),
                      deriv2(gf_gammaxy_, I, dx),
                      deriv2(gf_gammaxz_, I, dx),
                      deriv2(gf_gammayy_, I, dx),
                      deriv2(gf_gammayz_, I, dx),
                      deriv2(gf_gammazz_, I, dx),
                  },
                  //
                  deriv2(gf_alpha_, I, dx),
                  //
                  vec3<mat3<T, DN, DN>, UP>{
                      deriv2(gf_betax_, I, dx),
                      deriv2(gf_betay_, I, dx),
                      deriv2(gf_betaz_, I, dx),
                  })
  //
  {}
#endif
};

} // namespace Weyl

#endif // #ifndef WEYL_HXX
