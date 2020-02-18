#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include "tensor.hxx"

namespace Z4c {

template <typename T>
CCTK_ATTRIBUTE_ALWAYS_INLINE constexpr mat3<vec3<T, DN>, UP, UP>
calc_dgu(const mat3<T, UP, UP> &gu, const mat3<vec3<T, DN>, DN, DN> &dg) {
  // g_xy = g_xa g_yb g^ab
  // g_xy,c = (g_xa g_yb g^ab),c
  //        = g_xa,c g_yb g^ab + g_xa g_yb,c g^ab + g_xa g_yb g^ab,c
  // g_xa g_yb g^ab,c = g_xy,c - g_xa,c g_yb g^ab - g_xa g_yb,c g^ab
  //                  = g_xy,c - g_xy,c - g_xy,c
  //                  = - g_xy,c
  // g^ab,c = - g^ax g^by g_xy,c
  return mat3<vec3<T, DN>, UP, UP>(
      [&](int a, int b) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        return vec3<T, DN>([&](int c) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          return sum2(
              [&](int x, int y) { return -gu(a, x) * gu(b, y) * dg(x, y)(c); });
        });
      });
}

template <typename T>
CCTK_ATTRIBUTE_ALWAYS_INLINE constexpr mat3<vec3<T, DN>, UP, UP>
calc_dAu(const mat3<T, UP, UP> &gu, const mat3<vec3<T, DN>, UP, UP> &dgu,
         const mat3<T, DN, DN> &A, const mat3<vec3<T, DN>, DN, DN> &dA) {
  // A^ab,c = (g^ax g^by A_xy),c
  //        = g^ax,c g^by A_xy + g^ax g^by,c A_xy + g^ax g^by A_xy,c
  return mat3<vec3<T, DN>, UP, UP>(
      [&](int a, int b) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        return vec3<T, DN>([&](int c) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          return sum2([&](int x, int y) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return dgu(a, x)(c) * gu(b, y) * A(x, y)   //
                   + gu(a, x) * dgu(b, y)(c) * A(x, y) //
                   + gu(a, x) * gu(b, y) * dA(x, y)(c);
          });
        });
      });
}

template <typename T>
CCTK_ATTRIBUTE_ALWAYS_INLINE constexpr vec3<mat3<T, DN, DN>, DN>
calc_gammal(const mat3<vec3<T, DN>, DN, DN> dg) {
  // Gammal_abc
  return vec3<mat3<T, DN, DN>, DN>([&](int a) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return mat3<T, DN, DN>([&](int b, int c) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      return (dg(a, b)(c) + dg(a, c)(b) - dg(b, c)(a)) / 2;
    });
  });
}

template <typename T>
CCTK_ATTRIBUTE_ALWAYS_INLINE constexpr vec3<mat3<T, DN, DN>, UP>
calc_gamma(const mat3<T, UP, UP> &gu, const vec3<mat3<T, DN, DN>, DN> &Gammal) {
  // Gamma^a_bc
  return vec3<mat3<T, DN, DN>, UP>([&](int a) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    return mat3<T, DN, DN>([&](int b, int c) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      return sum1([&](int x) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        return gu(a, x) * Gammal(x)(b, c);
      });
    });
  });
}

} // namespace Z4c

#endif // #ifndef PHYSICS_HXX
