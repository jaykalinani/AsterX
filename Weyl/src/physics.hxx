#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include "tensor.hxx"

#include <dual.hxx>

#include <cmath>

namespace Weyl {
using namespace std;

template <typename T>
constexpr mat4<vec4<T, DN>, UP, UP>
calc_dgu(const mat4<T, UP, UP> &gu, const mat4<vec4<T, DN>, DN, DN> &dg) {
  // g_xy = g_xa g_yb g^ab
  // g_xy,c = (g_xa g_yb g^ab),c
  //        = g_xa,c g_yb g^ab + g_xa g_yb,c g^ab + g_xa g_yb g^ab,c
  // g_xa g_yb g^ab,c = g_xy,c - g_xa,c g_yb g^ab - g_xa g_yb,c g^ab
  //                  = g_xy,c - g_xy,c - g_xy,c
  //                  = - g_xy,c
  // g^ab,c = - g^ax g^by g_xy,c
  return mat4<vec4<T, DN>, UP, UP>([&](int a, int b) {
    return vec4<T, DN>([&](int c) {
      return sum42(
          [&](int x, int y) { return -gu(a, x) * gu(b, y) * dg(x, y)(c); });
    });
  });
}

template <typename T>
constexpr vec4<mat4<T, DN, DN>, DN>
calc_gammal(const mat4<vec4<T, DN>, DN, DN> dg) {
  // Gammal_abc
  return vec4<mat4<T, DN, DN>, DN>([&](int a) {
    return mat4<T, DN, DN>([&](int b, int c) {
      return (dg(a, b)(c) + dg(a, c)(b) - dg(b, c)(a)) / 2;
    });
  });
}

template <typename T>
constexpr vec4<mat4<T, DN, DN>, UP>
calc_gamma(const mat4<T, UP, UP> &gu, const vec4<mat4<T, DN, DN>, DN> &Gammal) {
  // Gamma^a_bc
  return vec4<mat4<T, DN, DN>, UP>([&](int a) {
    return mat4<T, DN, DN>([&](int b, int c) {
      return sum41([&](int x) { return gu(a, x) * Gammal(x)(b, c); });
    });
  });
}

template <typename T>
constexpr vec4<mat4<vec4<T, DN>, DN, DN>, DN>
calc_dgammal(const mat4<mat4<T, DN, DN>, DN, DN> &ddg) {
  // Gamma_abc,d
  return vec4<mat4<vec4<T, DN>, DN, DN>, DN>([&](int a) {
    return mat4<vec4<T, DN>, DN, DN>([&](int b, int c) {
      return vec4<T, DN>([&](int d) {
        return sum41([&](int x) {
          return (ddg(a, b)(c, d) + ddg(a, c)(b, d) - ddg(b, c)(a, d)) / 2;
        });
      });
    });
  });
}

template <typename T>
constexpr vec4<mat4<vec4<T, DN>, DN, DN>, UP>
calc_dgamma(const mat4<T, UP, UP> &gu, const mat4<vec4<T, DN>, UP, UP> &dgu,
            const vec4<mat4<T, DN, DN>, DN> &Gammal,
            const vec4<mat4<vec4<T, DN>, DN, DN>, DN> &dGammal) {
  // Gamma^a_bc,d
  return vec4<mat4<vec4<T, DN>, DN, DN>, UP>([&](int a) {
    return mat4<vec4<T, DN>, DN, DN>([&](int b, int c) {
      return vec4<T, DN>([&](int d) {
        return sum41([&](int x) {
          return dgu(a, x)(d) * Gammal(x)(b, c) +
                 gu(a, x) * dGammal(x)(b, c)(d);
        });
      });
    });
  });
}

template <typename T>
constexpr amat4<amat4<T, DN, DN>, DN, DN>
calc_riemann(const mat4<T, DN, DN> &g, const vec4<mat4<T, DN, DN>, UP> &Gamma,
             const vec4<mat4<vec4<T, DN>, DN, DN>, UP> &dGamma) {
  // Rm_abcd
  return amat4<amat4<T, DN, DN>, DN, DN>([&](int a, int b) {
    return amat4<T, DN, DN>([&](int c, int d) {
      return sum41([&](int x) {
        const T rmuxbcd =
            dGamma(x)(d, b)(c)                                              //
            - dGamma(x)(c, b)(d)                                            //
            + sum41([&](int y) { return Gamma(x)(c, y) * Gamma(y)(d, b); }) //
            - sum41([&](int y) { return Gamma(x)(d, y) * Gamma(y)(c, b); });
        return g(a, x) * rmuxbcd;
      });
    });
  });
}

template <typename T>
constexpr mat4<T, DN, DN>
calc_ricci(const mat4<T, UP, UP> &gu,
           const amat4<amat4<T, DN, DN>, DN, DN> &Rm) {
  // R_ab
  return mat4<T, DN, DN>([&](int a, int b) {
    return sum42([&](int x, int y) { return gu(x, y) * Rm(x, a)(y, b); });
  });
}

template <typename T>
constexpr amat4<amat4<T, DN, DN>, DN, DN>
calc_weyl(const mat4<T, DN, DN> &g, const amat4<amat4<T, DN, DN>, DN, DN> &Rm,
          const mat4<T, DN, DN> &R, const T &Rsc) {
  // C_abcd
  return amat4<amat4<T, DN, DN>, DN, DN>([&](int a, int b) {
    return amat4<T, DN, DN>([&](int c, int d) {
      return Rm(a, b)(c, d) //
             + 1 / T{4 - 2} *
                   (R(a, d) * g(b, c) - R(a, c) * g(b, d)    //
                    + R(b, c) * g(a, d) - R(b, d) * g(a, c)) //
             + 1 / T{(4 - 1) * (4 - 2)} * Rsc *
                   (g(a, c) * g(b, d) - g(a, d) * g(b, c));
    });
  });
}

template <typename T>
constexpr vec4<T, UP> raise(const mat4<T, UP, UP> &gu, const vec4<T, DN> &vl) {
  return vec4<T, UP>(
      [&](int a) { return sum41([&](int x) { return gu(a, x) * vl(x); }); });
}

template <typename T>
constexpr vec4<T, DN> lower(const mat4<T, DN, DN> &g, const vec4<T, UP> &v) {
  return vec4<T, DN>(
      [&](int a) { return sum41([&](int x) { return g(a, x) * v(x); }); });
}

template <typename T>
constexpr T dot(const vec4<T, DN> &vl, const vec4<T, UP> &v) {
  return sum41([&](int x) { return vl(x) * v(x); });
}

template <typename T>
constexpr vec4<T, UP> normalized(const mat4<T, DN, DN> &g,
                                 const vec4<T, UP> &v) {
  return v / sqrt(dot(lower(g, v), v));
}

template <typename T>
constexpr vec4<T, UP> projected(const mat4<T, DN, DN> &g, const vec4<T, UP> &v,
                                const vec4<T, UP> &w) {
  const auto z = zero<T>()();
  const auto wl = lower(g, w);
  const auto wlen2 = dot(wl, w);
  if (wlen2 < pow(1.0e-12, 2))
    return vec4<T, UP>([&](int a) { return z; });
  return dot(wl, v) / wlen2 * w;
}

template <typename T>
constexpr vec4<T, UP> projected1(const vec4<T, UP> &v, const vec4<T, DN> &wl,
                                 const vec4<T, UP> &w) {
  // assuming dot(wl, w) == 1
  return dot(wl, v) * w;
}

template <typename T>
constexpr vec4<T, UP> rejected(const mat4<T, DN, DN> &g, const vec4<T, UP> &v,
                               const vec4<T, UP> &w) {
  return v - projected(g, v, w);
}

template <typename T>
constexpr vec4<T, UP> rejected1(const vec4<T, UP> &v, const vec4<T, DN> &wl,
                                const vec4<T, UP> &w) {
  return v - projected1(v, wl, w);
}

template <typename T> constexpr vec4<T, UP> calc_et(const mat4<T, UP, UP> &gu) {
  const auto z = zero<T>()();
  const auto e = one<T>()();
  const vec4<T, DN> etl([&](int a) { return a == 0 ? e : z; });
  const auto et = raise(gu, etl);
  const auto etlen = sqrt(-dot(etl, et));
  assert(!isnan(etlen) && etlen >= 1.0e-12);
  return et / etlen;
}

template <typename T>
constexpr vec4<T, UP> calc_ephi(const vec4<T, UP> &x,
                                const mat4<T, DN, DN> &g) {
  const auto z = zero<T>()();
  vec4<T, UP> ephi(z, -x(2), x(1), z);
  const auto ephil = lower(g, ephi);
  const auto ephi_len2 = dot(ephil, ephi);
  if (ephi_len2 < pow(1.0e-12, 2))
    return vec4<T, UP>([&](int a) { return z; });
  ephi /= sqrt(ephi_len2);
  return ephi;
}

template <typename T>
constexpr vec4<T, UP> calc_etheta(const vec4<T, UP> &x,
                                  const mat4<T, DN, DN> &g,
                                  const vec4<T, UP> &ephi) {
  const auto z = zero<T>()();
  const T rho2 = pow(x(1), 2) + pow(x(2), 2);
  vec4<T, UP> etheta(z, x(1) * x(3), x(2) * x(3), -rho2);
  const auto ethetal = lower(g, etheta);
  const auto etheta_len2 = dot(ethetal, etheta);
  if (etheta_len2 < pow(1.0e-12, 2))
    return vec4<T, UP>([&](int a) { return z; });
  etheta /= sqrt(etheta_len2); // to improve accuracy
  etheta = rejected(g, etheta, ephi);
  etheta = normalized(g, etheta);
  return etheta;
}

template <typename T>
constexpr vec4<T, UP> calc_er(const vec4<T, UP> &x, const mat4<T, DN, DN> &g,
                              const vec4<T, UP> &etheta,
                              const vec4<T, UP> &ephi) {
  const auto z = zero<T>()();
  vec4<T, UP> er(z, x(1), x(2), x(3));
  const auto erl = lower(g, er);
  const auto er_len2 = dot(erl, er);
  if (er_len2 < pow(1.0e-12, 2))
    return vec4<T, UP>([&](int a) { return z; });
  er /= sqrt(er_len2); // to improve accuracy
  er = rejected(g, er, etheta);
  er = rejected(g, er, ephi);
  er = normalized(g, er);
  return er;
}

template <typename T>
constexpr vec4<vec4<T, DN>, UP>
calc_det(const mat4<T, UP, UP> &gu, const mat4<vec4<T, DN>, UP, UP> &dgu,
         const vec4<T, UP> &et, const vec4<mat4<T, DN, DN>, UP> &Gamma) {
  typedef vec4<T, DN> DT;
  typedef dual<T, DT> TDT;
  const mat4<TDT, UP, UP> gu1(
      [&](int a, int b) { return TDT(gu(a, b), dgu(a, b)); });
  const auto det1 = calc_et(gu1);
  const vec4<vec4<T, DN>, UP> det([&](int a) {
    return vec4<T, DN>([&](int b) {
      return det1(a).eps(b) +
             sum41([&](int x) { return Gamma(a)(b, x) * et(x); });
    });
  });
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b)
      assert(!isnan(det(a)(b)));
  return det;
}

template <typename T>
constexpr vec4<vec4<T, DN>, UP>
calc_dephi(const vec4<T, UP> &x, const mat4<T, DN, DN> &g,
           const mat4<vec4<T, DN>, DN, DN> &dg, const vec4<T, UP> &ephi,
           const vec4<mat4<T, DN, DN>, UP> &Gamma) {
  typedef vec4<T, DN> DT;
  typedef dual<T, DT> TDT;
  const vec4<TDT, UP> x1(
      [&](int a) { return TDT(x(a), DT([&](int b) { return T(a == b); })); });
  const mat4<TDT, DN, DN> g1(
      [&](int a, int b) { return TDT(g(a, b), dg(a, b)); });
  const auto dephi1 = calc_ephi(x1, g1);
  const vec4<vec4<T, DN>, UP> dephi([&](int a) {
    return vec4<T, DN>([&](int b) {
      return dephi1(a).eps(b) +
             sum41([&](int x) { return Gamma(a)(b, x) * ephi(x); });
    });
  });
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b)
      assert(!isnan(dephi(a)(b)));
  return dephi;
}

template <typename T>
constexpr vec4<vec4<T, DN>, UP>
calc_detheta(const vec4<T, UP> &x, const mat4<T, DN, DN> &g,
             const mat4<vec4<T, DN>, DN, DN> &dg, const vec4<T, UP> &ephi,
             const vec4<vec4<T, DN>, UP> &dephi, const vec4<T, UP> &etheta,
             const vec4<mat4<T, DN, DN>, UP> &Gamma) {
  typedef vec4<T, DN> DT;
  typedef dual<T, DT> TDT;
  const vec4<TDT, UP> x1(
      [&](int a) { return TDT(x(a), DT([&](int b) { return T(a == b); })); });
  const mat4<TDT, DN, DN> g1(
      [&](int a, int b) { return TDT(g(a, b), dg(a, b)); });
  const vec4<TDT, UP> ephi1([&](int a) { return TDT(ephi(a), dephi(a)); });
  const auto detheta1 = calc_etheta(x1, g1, ephi1);
  const vec4<vec4<T, DN>, UP> detheta([&](int a) {
    return vec4<T, DN>([&](int b) {
      return detheta1(a).eps(b) +
             sum41([&](int x) { return Gamma(a)(b, x) * etheta(x); });
    });
  });
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b)
      assert(!isnan(detheta(a)(b)));
  return detheta;
}

template <typename T>
constexpr vec4<vec4<T, DN>, UP>
calc_der(const vec4<T, UP> &x, const mat4<T, DN, DN> &g,
         const mat4<vec4<T, DN>, DN, DN> &dg, const vec4<T, UP> &ephi,
         const vec4<vec4<T, DN>, UP> &dephi, const vec4<T, UP> &etheta,
         const vec4<vec4<T, DN>, UP> &detheta, const vec4<T, UP> &er,
         const vec4<mat4<T, DN, DN>, UP> &Gamma) {
  typedef vec4<T, DN> DT;
  typedef dual<T, DT> TDT;
  const vec4<TDT, UP> x1(
      [&](int a) { return TDT(x(a), DT([&](int b) { return T(a == b); })); });
  const mat4<TDT, DN, DN> g1(
      [&](int a, int b) { return TDT(g(a, b), dg(a, b)); });
  const vec4<TDT, UP> ephi1([&](int a) { return TDT(ephi(a), dephi(a)); });
  const vec4<TDT, UP> etheta1(
      [&](int a) { return TDT(etheta(a), detheta(a)); });
  const auto der1 = calc_er(x1, g1, ephi1, etheta1);
  const vec4<vec4<T, DN>, UP> der([&](int a) {
    return vec4<T, DN>([&](int b) {
      return der1(a).eps(b) +
             sum41([&](int x) { return Gamma(a)(b, x) * er(x); });
    });
  });
  for (int a = 0; a < 4; ++a)
    for (int b = 0; b < 4; ++b)
      assert(!isnan(der(a)(b)));
  return der;
}

} // namespace Weyl

#endif // #ifndef PHYSICS_HXX
