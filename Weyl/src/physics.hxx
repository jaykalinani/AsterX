#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include "tensor.hxx"

namespace Weyl {

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

#if 0
template <typename T>
 constexpr mat4<vec4<T, DN>, UP, UP>
calc_dAu(const mat4<T, UP, UP> &gu, const mat4<vec4<T, DN>, UP, UP> &dgu,
         const mat4<T, DN, DN> &A, const mat4<vec4<T, DN>, DN, DN> &dA) {
  // A^ab,c = (g^ax g^by A_xy),c
  //        = g^ax,c g^by A_xy + g^ax g^by,c A_xy + g^ax g^by A_xy,c
  return mat4<vec4<T, DN>, UP, UP>(
      [&](int a, int b)  {
        return vec4<T, DN>([&](int c)  {
          return sum42([&](int x, int y)  {
            return dgu(a, x)(c) * gu(b, y) * A(x, y)   //
                   + gu(a, x) * dgu(b, y)(c) * A(x, y) //
                   + gu(a, x) * gu(b, y) * dA(x, y)(c);
          });
        });
      });
}
#endif

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
calc_riemann(const vec4<mat4<T, DN, DN>, DN> &Gammal,
             const vec4<mat4<T, DN, DN>, UP> &Gamma,
             const vec4<mat4<vec4<T, DN>, DN, DN>, DN> &dGammal) {
  // Rm_abcd
  return amat4<amat4<T, DN, DN>, DN, DN>([&](int a, int b) {
    return amat4<T, DN, DN>([&](int c, int d) {
      return dGammal(a)(d, b)(c)                                              //
             - dGammal(a)(c, b)(d)                                            //
             + sum41([&](int x) { return Gammal(a)(c, x) * Gamma(x)(d, b); }) //
             - sum41([&](int x) { return Gammal(a)(d, x) * Gamma(x)(c, b); });
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
      return Rm(a, b)(c, d) +
             1 / T{4 - 2} *
                 (R(a, d) * g(b, c) - R(a, c) * g(b, d) + R(b, c) * g(a, d) -
                  R(b, d) * g(a, c)) +
             1 / T{(4 - 1) * (4 - 2)} * Rsc *
                 (g(a, c) * g(b, d) - g(a, d) * g(b, d));
    });
  });
}

} // namespace Weyl

#endif // #ifndef PHYSICS_HXX
