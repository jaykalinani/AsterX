#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include "tensor.hxx"

namespace Z4c {

template <typename T>
Z4C_INLINE constexpr mat3<vec3<T, DN>, UP, UP>
calc_dgu(const mat3<T, UP, UP> &gu, const mat3<vec3<T, DN>, DN, DN> &dg) {
  // g_xy = g_xa g_yb g^ab
  // g_xy,c = (g_xa g_yb g^ab),c
  //        = g_xa,c g_yb g^ab + g_xa g_yb,c g^ab + g_xa g_yb g^ab,c
  // g_xa g_yb g^ab,c = g_xy,c - g_xa,c g_yb g^ab - g_xa g_yb,c g^ab
  //                  = g_xy,c - g_xy,c - g_xy,c
  //                  = - g_xy,c
  // g^ab,c = - g^ax g^by g_xy,c
  return mat3<vec3<T, DN>, UP, UP>([&] Z4C_INLINE(int a, int b) {
    return vec3<T, DN>([&] Z4C_INLINE(int c) {
      // return sum2sym([&]Z4C_INLINE(int x, int y) {
      //   return -gu(a, x) * gu(b, y) * dg(x, y)(c);
      // });
      return sum1([&] Z4C_INLINE(int x) {
        return -gu(a, x) *
               sum1([&] Z4C_INLINE(int y) { return gu(b, y) * dg(x, y)(c); });
      });
    });
  });
}

template <typename T>
Z4C_INLINE constexpr mat3<vec3<T, DN>, UP, UP>
calc_dAu(const mat3<T, UP, UP> &gu, const mat3<vec3<T, DN>, UP, UP> &dgu,
         const mat3<T, DN, DN> &A, const mat3<vec3<T, DN>, DN, DN> &dA) {
  // A^ab,c = (g^ax g^by A_xy),c
  //        = g^ax,c g^by A_xy + g^ax g^by,c A_xy + g^ax g^by A_xy,c
  return mat3<vec3<T, DN>, UP, UP>([&] Z4C_INLINE(int a, int b) {
    return vec3<T, DN>([&] Z4C_INLINE(int c) {
      // return sum2sym([&]Z4C_INLINE(int x, int y) {
      //   return dgu(a, x)(c) * gu(b, y) * A(x, y)   //
      //          + gu(a, x) * dgu(b, y)(c) * A(x, y) //
      //          + gu(a, x) * gu(b, y) * dA(x, y)(c);
      // });
      return sum1([&] Z4C_INLINE(int x) {
        return gu(b, x) * sum1([&] Z4C_INLINE(int y) {
                 return dgu(a, y)(c) * A(x, y)   //
                        + dgu(b, y)(c) * A(x, y) //
                        + gu(b, y) * dA(x, y)(c);
               });
      });
    });
  });
}

template <typename T>
Z4C_INLINE constexpr vec3<mat3<T, DN, DN>, DN>
calc_gammal(const mat3<vec3<T, DN>, DN, DN> &dg) {
  // Gammal_abc
  return vec3<mat3<T, DN, DN>, DN>([&] Z4C_INLINE(int a) {
    return mat3<T, DN, DN>([&] Z4C_INLINE(int b, int c) {
      return (dg(a, b)(c) + dg(a, c)(b) - dg(b, c)(a)) / 2;
    });
  });
}

template <typename T>
Z4C_INLINE constexpr vec3<mat3<T, DN, DN>, UP>
calc_gamma(const mat3<T, UP, UP> &gu, const vec3<mat3<T, DN, DN>, DN> &Gammal) {
  // Gamma^a_bc
  return vec3<mat3<T, DN, DN>, UP>([&] Z4C_INLINE(int a) {
    return mat3<T, DN, DN>([&] Z4C_INLINE(int b, int c) {
      return sum1([&] Z4C_INLINE(int x) { return gu(a, x) * Gammal(x)(b, c); });
    });
  });
}

template <typename T>
Z4C_INLINE constexpr vec3<vec3<vec3<T, UP>, DN>, DN>
calc_gammalu(const mat3<T, UP, UP> &gu,
             const vec3<mat3<T, DN, DN>, DN> &Gammal) {
  // Gamma_ab^c
  return vec3<vec3<vec3<T, UP>, DN>, DN>([&] Z4C_INLINE(int a) {
    return vec3<vec3<T, UP>, DN>([&] Z4C_INLINE(int b) {
      return vec3<T, UP>([&] Z4C_INLINE(int c) {
        return sum1(
            [&] Z4C_INLINE(int x) { return Gammal(a)(b, x) * gu(x, c); });
      });
    });
  });
}

} // namespace Z4c

#endif // #ifndef PHYSICS_HXX
