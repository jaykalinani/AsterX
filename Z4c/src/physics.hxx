#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include "tensor.hxx"

namespace Z4c {

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<vec3<T> >
calc_dgu(const mat3<T> &gu, const mat3<vec3<T> > &dg) {
  // g_xy = g_xa g_yb g^ab
  // g_xy,c = (g_xa g_yb g^ab),c
  //        = g_xa,c g_yb g^ab + g_xa g_yb,c g^ab + g_xa g_yb g^ab,c
  // g_xa g_yb g^ab,c = g_xy,c - g_xa,c g_yb g^ab - g_xa g_yb,c g^ab
  //                  = g_xy,c - g_xy,c - g_xy,c
  //                  = - g_xy,c
  // g^ab,c = - g^ax g^by g_xy,c
  return mat3<vec3<T> >([&](int a, int b) {
    return vec3<T>([&](int c) {
      return sum2(
          [&](int x, int y) { return -gu(a, x) * gu(b, y) * dg(x, y)(c); });
    });
  });
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<vec3<T> >
calc_dAu(const mat3<T> &gu, const mat3<vec3<T> > &dgu, const mat3<T> &A,
         const mat3<vec3<T> > &dA) {
  // A^ab,c = (g^ax g^by A_xy),c
  //        = g^ax,c g^by A_xy + g^ax g^by,c A_xy + g^ax g^by A_xy,c
  return mat3<vec3<T> >([&](int a, int b) {
    return vec3<T>([&](int c) {
      return sum2([&](int x, int y) {
        return dgu(a, x)(c) * gu(b, y) * A(x, y)   //
               + gu(a, x) * dgu(b, y)(c) * A(x, y) //
               + gu(a, x) * gu(b, y) * dA(x, y)(c);
      });
    });
  });
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3<mat3<T> >
calc_gammal(const mat3<vec3<T> > dg) {
  // Gammal_abc
  return vec3<mat3<T> >([&](int a) {
    return mat3<T>([&](int b, int c) {
      return (dg(a, b)(c) + dg(a, c)(b) - dg(b, c)(a)) / 2;
    });
  });
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3<mat3<T> >
calc_gamma(const mat3<T> &gu, const vec3<mat3<T> > &Gammal) {
  // Gamma^a_bc
  return vec3<mat3<T> >([&](int a) {
    return mat3<T>([&](int b, int c) {
      return sum1([&](int x) { return gu(a, x) * Gammal(x)(b, c); });
    });
  });
}

} // namespace Z4c

#endif // #ifndef PHYSICS_HXX
