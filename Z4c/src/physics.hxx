#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include "tensor.hxx"

namespace Z4c {

template <typename T>
mat3<vec3<T> > calc_dgu(const mat3<T> &gu, const mat3<vec3<T> > &dg) {
  // g_xy = g_xa g_yb g^ab
  // g_xy,c = (g_xa g_yb g^ab),c
  //        = g_xa,c g_yb g^ab + g_xa g_yb,c g^ab + g_xa g_yb g^ab,c
  // g_xa g_yb g^ab,c = g_xy,c - g_xa,c g_yb g^ab - g_xa g_yb,c g^ab
  //                  = g_xy,c - g_xy,c - g_xy,c
  //                  = - g_xy,c
  // g^ab,c = - g^ax g^by g_xy,c
  mat3<vec3<T> > dgu;
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      for (int c = 0; c < 3; ++c) {
        T s = 0;
        for (int x = 0; x < 3; ++x)
          for (int y = 0; y < 3; ++y)
            s -= gu(a, x) * gu(b, y) * dg(x, y)(c);
        dgu(a, b)(c) = s;
      }
  return dgu;
}

template <typename T>
mat3<vec3<T> > calc_dAu(const mat3<T> &gu, const mat3<vec3<T> > &dgu,
                        const mat3<T> &A, const mat3<vec3<T> > &dA) {
  // A^ab,c = (g^ax g^by A_xy),c
  //        = g^ax,c g^by A_xy + g^ax g^by,c A_xy + g^ax g^by A_xy,c
  mat3<vec3<T> > dAu;
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      for (int c = 0; c < 3; ++c) {
        T s = 0;
        for (int x = 0; x < 3; ++x)
          for (int y = 0; y < 3; ++y)
            s += dgu(a, x)(c) * gu(b, y) * A(x, y) +
                 gu(a, x) * dgu(b, y)(c) * A(x, y) +
                 gu(a, x) * gu(b, y) * dA(x, y)(c);
        dAu(a, b)(c) = s;
      }
  return dAu;
}

template <typename T>
void calc_gamma(const mat3<T> &g, const mat3<T> &gu, const mat3<vec3<T> > dg,
                vec3<mat3<T> > &restrict Gammal, vec3<mat3<T> > &restrict Gamma,
                vec3<T> &restrict Gam) {
  // Gammal_abc
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      for (int c = b; c < 3; ++c)
        Gammal(a)(b, c) = (dg(a, b)(c) + dg(a, c)(b) - dg(a, b)(c)) / 2;

  // Gamma^a_bc
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      for (int c = b; c < 3; ++c) {
        T s = 0;
        for (int x = 0; x < 3; ++x)
          s += gu(a, x) * Gammal(x)(b, c);
        Gamma(a)(b, c) = s;
      }

  // Gam^a
  for (int a = 0; a < 3; ++a) {
    T s = 0;
    for (int x = 0; x < 3; ++x)
      for (int y = 0; y < 3; ++y)
        s += gu(x, y) * Gamma(a)(x, y);
    Gam(a) = s;
  }
}

template <typename T>
mat3<T> calc_ricci(const T &chi, const vec3<T> &dchi, const mat3<T> &DDchi,
                   const mat3<T> &gammat, const mat3<T> &gammatu,
                   const mat3<mat3<T> > &ddgammat,
                   const vec3<mat3<T> > &Gammatdl,
                   const vec3<mat3<T> > &Gammatd, const vec3<T> &Gamtd,
                   const vec3<vec3<T> > &dGamt) {
  // (8)
  mat3<T> Rchi;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) {
      T s = 0;
      s += DDchi(a, b) / (2 * chi);
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
          s += gammat(a, b) * gammatu(x, y) * DDchi(x, y);
      s -= dchi(a) * dchi(b) / (4 * chi * chi);
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
          s -= gammat(a, b) * gammatu(x, y) * dchi(x) * dchi(y) * 3 /
               (4 * chi * chi);
      Rchi(a, b) = s;
    }

  // (9)
  mat3<T> Rt;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b) {
      T s = 0;
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
          s -= gammatu(x, y) * ddgammat(a, b)(x, y) / 2;
      for (int x = 0; x < 3; ++x)
        s += gammat(x, a) * dGamt(x)(b) + gammat(x, b) * dGamt(x)(a);
      for (int x = 0; x < 3; ++x)
        s += Gamtd(x) * Gammatdl(a)(b, x) + Gamtd(x) * Gammatdl(b)(a, x);
      for (int x = 0; x < 3; ++x)
        for (int y = 0; y < 3; ++y)
          for (int z = 0; z < 3; ++z)
            s += gammatu(y, z) * (2 * Gammatd(x)(y, a) * Gammatdl(b)(x, z) +
                                  2 * Gammatd(x)(y, b) * Gammatdl(a)(x, z) +
                                  Gammatd(x)(a, z) * Gammatdl(x)(y, b) +
                                  Gammatd(x)(b, z) * Gammatdl(x)(y, a));
      Rt(a, b) = s;
    }

  // (7)
  mat3<T> R;
  for (int a = 0; a < 3; ++a)
    for (int b = 0; b < 3; ++b)
      R(a, b) = Rchi(a, b) + Rt(a, b);

  return R;
}

} // namespace Z4c

#endif // #ifndef PHYSICS_HXX
