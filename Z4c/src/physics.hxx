#ifndef PHYSICS_HXX
#define PHYSICS_HXX

#include <defs.hxx>
#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>

namespace Z4c {
using namespace Arith;

#if 0
template <typename T, int D, dnup_t dnup1, dnup_t dnup2, symm_t symm>
gmat<T, D - 1, dnup1, dnup2, symm>
calc_submatrix(const gmat<T, D, dnup1, dnup2, symm> A, int i0, int j0) {
  return gmat<T, D - 1, dnup1, dnup2, symm>([&](int is, int js) {
    const int i = is + (is >= i0);
    const int j = js + (js >= j0);
    return A(i, j);
  });
}

template <typename T, int D, dnup_t dnup1, dnup_t dnup2, symm_t symm>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const gmat<T, D, dnup1, dnup2, symm> &g) {
  if constexpr (D == 0) {
    return one<T>();
  } else {
    T detg = zero<T>();
    for (int i = 0; i < D; ++i)
      detg += bitsign(i) * calc_det(calc_submatrix(g, i, 0));
    return detg;
  }
}

template <typename T, int D, dnup_t dnup1, dnup_t dnup2, symm_t symm>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr gmat<T, D, !dnup1, !dnup2, symm>
calc_inv(const gmat<T, D, dnup1, dnup2, symm> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return gmat<T, D, !dnup1, !dnup2, symm>([&](int i, int j) {
    return bitsign(i) * bitsign(j) * detg1 * calc_det(calc_submatrix(g, j, i));
  });
}
#endif

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const smat<T, 0, dnup1, dnup2> &g) {
  return one<T>();
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const smat<T, 1, dnup1, dnup2> &g) {
  return g(0, 0);
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const smat<T, 2, dnup1, dnup2> &g) {
  return g(0, 0) * g(1, 1) - pow2(g(0, 1));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const smat<T, 3, dnup1, dnup2> &g) {
  return 2 * g(0, 1) * g(0, 2) * g(1, 2) - g(2, 2) * pow2(g(0, 1)) -
         g(1, 1) * pow2(g(0, 2)) +
         g(0, 0) * (g(1, 1) * g(2, 2) - pow2(g(1, 2)));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const smat<T, 4, dnup1, dnup2> &g) {
  return 2 * g(0, 0) * g(1, 2) * g(1, 3) * g(2, 3) -
         2 * g(0, 3) *
             (g(0, 2) * (g(1, 2) * g(1, 3) - g(1, 1) * g(2, 3)) +
              g(0, 1) * (-(g(1, 3) * g(2, 2)) + g(1, 2) * g(2, 3))) +
         g(0, 0) * g(1, 1) * g(2, 2) * g(3, 3) +
         2 * g(0, 1) * g(0, 2) * (-(g(1, 3) * g(2, 3)) + g(1, 2) * g(3, 3)) -
         g(2, 2) * g(3, 3) * pow2(g(0, 1)) - g(0, 0) * g(3, 3) * pow2(g(1, 2)) +
         pow2(g(0, 3)) * (-(g(1, 1) * g(2, 2)) + pow2(g(1, 2))) -
         g(0, 0) * g(2, 2) * pow2(g(1, 3)) +
         pow2(g(0, 2)) * (-(g(1, 1) * g(3, 3)) + pow2(g(1, 3))) -
         g(0, 0) * g(1, 1) * pow2(g(2, 3)) + pow2(g(0, 1)) * pow2(g(2, 3));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const mat<T, 0, dnup1, dnup2> &g) {
  return one<T>();
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const mat<T, 1, dnup1, dnup2> &g) {
  return g(0, 0);
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const mat<T, 2, dnup1, dnup2> &g) {
  return -(g(0, 1) * g(1, 0)) + g(0, 0) * g(1, 1);
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const mat<T, 3, dnup1, dnup2> &g) {
  return g(0, 2) * (-(g(1, 1) * g(2, 0)) + g(1, 0) * g(2, 1)) +
         g(0, 1) * (g(1, 2) * g(2, 0) - g(1, 0) * g(2, 2)) +
         g(0, 0) * (-(g(1, 2) * g(2, 1)) + g(1, 1) * g(2, 2));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T
calc_det(const mat<T, 4, dnup1, dnup2> &g) {
  return g(0, 1) * g(1, 3) * g(2, 2) * g(3, 0) -
         g(0, 1) * g(1, 2) * g(2, 3) * g(3, 0) -
         g(0, 0) * g(1, 3) * g(2, 2) * g(3, 1) +
         g(0, 0) * g(1, 2) * g(2, 3) * g(3, 1) -
         g(0, 1) * g(1, 3) * g(2, 0) * g(3, 2) +
         g(0, 0) * g(1, 3) * g(2, 1) * g(3, 2) +
         g(0, 1) * g(1, 0) * g(2, 3) * g(3, 2) -
         g(0, 0) * g(1, 1) * g(2, 3) * g(3, 2) +
         g(0, 3) * (g(1, 2) * (g(2, 1) * g(3, 0) - g(2, 0) * g(3, 1)) +
                    g(1, 1) * (-(g(2, 2) * g(3, 0)) + g(2, 0) * g(3, 2)) +
                    g(1, 0) * (g(2, 2) * g(3, 1) - g(2, 1) * g(3, 2))) +
         g(0, 1) * g(1, 2) * g(2, 0) * g(3, 3) -
         g(0, 0) * g(1, 2) * g(2, 1) * g(3, 3) -
         g(0, 1) * g(1, 0) * g(2, 2) * g(3, 3) +
         g(0, 0) * g(1, 1) * g(2, 2) * g(3, 3) +
         g(0, 2) * (g(1, 3) * (-(g(2, 1) * g(3, 0)) + g(2, 0) * g(3, 1)) +
                    g(1, 1) * (g(2, 3) * g(3, 0) - g(2, 0) * g(3, 3)) +
                    g(1, 0) * (-(g(2, 3) * g(3, 1)) + g(2, 1) * g(3, 3)));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 0, !dnup1, !dnup2>
calc_inv(const smat<T, 0, dnup1, dnup2> &g, const T &detg) {
  return smat<T, 0, !dnup1, !dnup2>([&](int i, int j) { return one<T>()(); });
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 1, !dnup1, !dnup2>
calc_inv(const smat<T, 1, dnup1, dnup2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return smat<T, 1, !dnup1, !dnup2>{detg1};
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 2, !dnup1, !dnup2>
calc_inv(const smat<T, 2, dnup1, dnup2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * smat<T, 2, !dnup1, !dnup2>{
                     g(1, 1),
                     -g(0, 1),

                     g(0, 0),
                 };
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 3, !dnup1, !dnup2>
calc_inv(const smat<T, 3, dnup1, dnup2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * smat<T, 3, !dnup1, !dnup2>{
                     g(1, 1) * g(2, 2) - pow2(g(1, 2)),
                     g(0, 2) * g(1, 2) - g(0, 1) * g(2, 2),
                     -(g(0, 2) * g(1, 1)) + g(0, 1) * g(1, 2),
                     g(0, 0) * g(2, 2) - pow2(g(0, 2)),
                     g(0, 1) * g(0, 2) - g(0, 0) * g(1, 2),
                     g(0, 0) * g(1, 1) - pow2(g(0, 1)),
                 };
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<T, 4, !dnup1, !dnup2>
calc_inv(const smat<T, 4, dnup1, dnup2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * smat<T, 4, !dnup1, !dnup2>{
                     2 * g(1, 2) * g(1, 3) * g(2, 3) - g(3, 3) * pow2(g(1, 2)) -
                         g(2, 2) * pow2(g(1, 3)) +
                         g(1, 1) * (g(2, 2) * g(3, 3) - pow2(g(2, 3))),
                     g(0, 3) * (g(1, 3) * g(2, 2) - g(1, 2) * g(2, 3)) +
                         g(0, 2) * (-(g(1, 3) * g(2, 3)) + g(1, 2) * g(3, 3)) +
                         g(0, 1) * (-(g(2, 2) * g(3, 3)) + pow2(g(2, 3))),
                     g(0, 3) * (-(g(1, 2) * g(1, 3)) + g(1, 1) * g(2, 3)) +
                         g(0, 1) * (-(g(1, 3) * g(2, 3)) + g(1, 2) * g(3, 3)) +
                         g(0, 2) * (-(g(1, 1) * g(3, 3)) + pow2(g(1, 3))),
                     g(0, 2) * (-(g(1, 2) * g(1, 3)) + g(1, 1) * g(2, 3)) +
                         g(0, 1) * (g(1, 3) * g(2, 2) - g(1, 2) * g(2, 3)) +
                         g(0, 3) * (-(g(1, 1) * g(2, 2)) + pow2(g(1, 2))),

                     2 * g(0, 2) * g(0, 3) * g(2, 3) - g(3, 3) * pow2(g(0, 2)) -
                         g(2, 2) * pow2(g(0, 3)) +
                         g(0, 0) * (g(2, 2) * g(3, 3) - pow2(g(2, 3))),
                     -(g(0, 3) * (g(0, 2) * g(1, 3) + g(0, 1) * g(2, 3))) +
                         g(0, 1) * g(0, 2) * g(3, 3) +
                         g(0, 0) * (g(1, 3) * g(2, 3) - g(1, 2) * g(3, 3)) +
                         g(1, 2) * pow2(g(0, 3)),
                     g(0, 1) * g(0, 3) * g(2, 2) -
                         g(0, 2) * (g(0, 3) * g(1, 2) + g(0, 1) * g(2, 3)) +
                         g(0, 0) * (-(g(1, 3) * g(2, 2)) + g(1, 2) * g(2, 3)) +
                         g(1, 3) * pow2(g(0, 2)),

                     2 * g(0, 1) * g(0, 3) * g(1, 3) - g(3, 3) * pow2(g(0, 1)) -
                         g(1, 1) * pow2(g(0, 3)) +
                         g(0, 0) * (g(1, 1) * g(3, 3) - pow2(g(1, 3))),
                     -(g(0, 1) * g(0, 3) * g(1, 2)) +
                         g(0, 2) * (g(0, 3) * g(1, 1) - g(0, 1) * g(1, 3)) +
                         g(0, 0) * (g(1, 2) * g(1, 3) - g(1, 1) * g(2, 3)) +
                         g(2, 3) * pow2(g(0, 1)),

                     2 * g(0, 1) * g(0, 2) * g(1, 2) - g(2, 2) * pow2(g(0, 1)) -
                         g(1, 1) * pow2(g(0, 2)) +
                         g(0, 0) * (g(1, 1) * g(2, 2) - pow2(g(1, 2))),
                 };
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 0, !dnup1, !dnup2>
calc_inv(const mat<T, 0, dnup1, dnup2> &g, const T &detg) {
  return mat<T, 0, !dnup1, !dnup2>([&](int i, int j) { return one<T>()(); });
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 1, !dnup1, !dnup2>
calc_inv(const mat<T, 1, dnup1, dnup2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return mat<T, 1, !dnup1, !dnup2>{detg1};
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 2, !dnup1, !dnup2>
calc_inv(const mat<T, 2, dnup1, dnup2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * mat<T, 2, !dnup1, !dnup2>{
                     g(1, 1),
                     -g(0, 1),

                     -g(1, 0),
                     g(0, 0),
                 };
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 3, !dnup1, !dnup2>
calc_inv(const mat<T, 3, dnup1, dnup2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * mat<T, 3, !dnup1, !dnup2>{
                     -(g(1, 2) * g(2, 1)) + g(1, 1) * g(2, 2),
                     g(0, 2) * g(2, 1) - g(0, 1) * g(2, 2),
                     -(g(0, 2) * g(1, 1)) + g(0, 1) * g(1, 2),

                     g(1, 2) * g(2, 0) - g(1, 0) * g(2, 2),
                     -(g(0, 2) * g(2, 0)) + g(0, 0) * g(2, 2),
                     g(0, 2) * g(1, 0) - g(0, 0) * g(1, 2),

                     -(g(1, 1) * g(2, 0)) + g(1, 0) * g(2, 1),
                     g(0, 1) * g(2, 0) - g(0, 0) * g(2, 1),
                     -(g(0, 1) * g(1, 0)) + g(0, 0) * g(1, 1),
                 };
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr mat<T, 4, !dnup1, !dnup2>
calc_inv(const mat<T, 4, dnup1, dnup2> &g, const T &detg) {
  const T detg1 = 1 / detg;
  return detg1 * mat<T, 4, !dnup1, !dnup2>{
                     g(1, 3) * (-(g(2, 2) * g(3, 1)) + g(2, 1) * g(3, 2)) +
                         g(1, 2) * (g(2, 3) * g(3, 1) - g(2, 1) * g(3, 3)) +
                         g(1, 1) * (-(g(2, 3) * g(3, 2)) + g(2, 2) * g(3, 3)),
                     g(0, 3) * (g(2, 2) * g(3, 1) - g(2, 1) * g(3, 2)) +
                         g(0, 2) * (-(g(2, 3) * g(3, 1)) + g(2, 1) * g(3, 3)) +
                         g(0, 1) * (g(2, 3) * g(3, 2) - g(2, 2) * g(3, 3)),
                     g(0, 3) * (-(g(1, 2) * g(3, 1)) + g(1, 1) * g(3, 2)) +
                         g(0, 2) * (g(1, 3) * g(3, 1) - g(1, 1) * g(3, 3)) +
                         g(0, 1) * (-(g(1, 3) * g(3, 2)) + g(1, 2) * g(3, 3)),
                     g(0, 3) * (g(1, 2) * g(2, 1) - g(1, 1) * g(2, 2)) +
                         g(0, 2) * (-(g(1, 3) * g(2, 1)) + g(1, 1) * g(2, 3)) +
                         g(0, 1) * (g(1, 3) * g(2, 2) - g(1, 2) * g(2, 3)),

                     g(1, 3) * (g(2, 2) * g(3, 0) - g(2, 0) * g(3, 2)) +
                         g(1, 2) * (-(g(2, 3) * g(3, 0)) + g(2, 0) * g(3, 3)) +
                         g(1, 0) * (g(2, 3) * g(3, 2) - g(2, 2) * g(3, 3)),
                     g(0, 3) * (-(g(2, 2) * g(3, 0)) + g(2, 0) * g(3, 2)) +
                         g(0, 2) * (g(2, 3) * g(3, 0) - g(2, 0) * g(3, 3)) +
                         g(0, 0) * (-(g(2, 3) * g(3, 2)) + g(2, 2) * g(3, 3)),
                     g(0, 3) * (g(1, 2) * g(3, 0) - g(1, 0) * g(3, 2)) +
                         g(0, 2) * (-(g(1, 3) * g(3, 0)) + g(1, 0) * g(3, 3)) +
                         g(0, 0) * (g(1, 3) * g(3, 2) - g(1, 2) * g(3, 3)),
                     g(0, 3) * (-(g(1, 2) * g(2, 0)) + g(1, 0) * g(2, 2)) +
                         g(0, 2) * (g(1, 3) * g(2, 0) - g(1, 0) * g(2, 3)) +
                         g(0, 0) * (-(g(1, 3) * g(2, 2)) + g(1, 2) * g(2, 3)),

                     g(1, 3) * (-(g(2, 1) * g(3, 0)) + g(2, 0) * g(3, 1)) +
                         g(1, 1) * (g(2, 3) * g(3, 0) - g(2, 0) * g(3, 3)) +
                         g(1, 0) * (-(g(2, 3) * g(3, 1)) + g(2, 1) * g(3, 3)),
                     g(0, 3) * (g(2, 1) * g(3, 0) - g(2, 0) * g(3, 1)) +
                         g(0, 1) * (-(g(2, 3) * g(3, 0)) + g(2, 0) * g(3, 3)) +
                         g(0, 0) * (g(2, 3) * g(3, 1) - g(2, 1) * g(3, 3)),
                     g(0, 3) * (-(g(1, 1) * g(3, 0)) + g(1, 0) * g(3, 1)) +
                         g(0, 1) * (g(1, 3) * g(3, 0) - g(1, 0) * g(3, 3)) +
                         g(0, 0) * (-(g(1, 3) * g(3, 1)) + g(1, 1) * g(3, 3)),
                     g(0, 3) * (g(1, 1) * g(2, 0) - g(1, 0) * g(2, 1)) +
                         g(0, 1) * (-(g(1, 3) * g(2, 0)) + g(1, 0) * g(2, 3)) +
                         g(0, 0) * (g(1, 3) * g(2, 1) - g(1, 1) * g(2, 3)),

                     g(1, 2) * (g(2, 1) * g(3, 0) - g(2, 0) * g(3, 1)) +
                         g(1, 1) * (-(g(2, 2) * g(3, 0)) + g(2, 0) * g(3, 2)) +
                         g(1, 0) * (g(2, 2) * g(3, 1) - g(2, 1) * g(3, 2)),
                     g(0, 2) * (-(g(2, 1) * g(3, 0)) + g(2, 0) * g(3, 1)) +
                         g(0, 1) * (g(2, 2) * g(3, 0) - g(2, 0) * g(3, 2)) +
                         g(0, 0) * (-(g(2, 2) * g(3, 1)) + g(2, 1) * g(3, 2)),
                     g(0, 2) * (g(1, 1) * g(3, 0) - g(1, 0) * g(3, 1)) +
                         g(0, 1) * (-(g(1, 2) * g(3, 0)) + g(1, 0) * g(3, 2)) +
                         g(0, 0) * (g(1, 2) * g(3, 1) - g(1, 1) * g(3, 2)),
                     g(0, 2) * (-(g(1, 1) * g(2, 0)) + g(1, 0) * g(2, 1)) +
                         g(0, 1) * (g(1, 2) * g(2, 0) - g(1, 0) * g(2, 2)) +
                         g(0, 0) * (-(g(1, 2) * g(2, 1)) + g(1, 1) * g(2, 2)),
                 };
}

template <typename T, int D, dnup_t dnup1, dnup_t dnup2>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr T calc_trace(
    const smat<T, D, dnup1, dnup2> &A, const smat<T, D, !dnup1, !dnup2> &gu) {
  return sum_symm<D>([&](int x, int y)
                         ARITH_INLINE { return gu(x, y) * A(x, y); });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<vec<T, D, DN>, D, UP, UP>
calc_dgu(const smat<T, D, UP, UP> &gu,
         const smat<vec<T, D, DN>, D, DN, DN> &dg) {
  // g_xy = g_xa g_yb g^ab
  // g_xy,c = (g_xa g_yb g^ab),c
  //        = g_xa,c g_yb g^ab + g_xa g_yb,c g^ab + g_xa g_yb g^ab,c
  // g_xa g_yb g^ab,c = g_xy,c - g_xa,c g_yb g^ab - g_xa g_yb,c g^ab
  //                  = g_xy,c - g_xy,c - g_xy,c
  //                  = - g_xy,c
  // g^ab,c = - g^ax g^by g_xy,c
  return smat<vec<T, D, DN>, D, UP, UP>([&](int a, int b) ARITH_INLINE {
    return vec<T, D, DN>([&](int c) ARITH_INLINE {
      // return sum2sym([&](int x, int y) ARITH_INLINE {
      //   return -gu(a, x) * gu(b, y) * dg(x, y)(c);
      // });
      return sum<D>([&](int x) ARITH_INLINE {
        return -gu(a, x) * sum<D>([&](int y) ARITH_INLINE {
          return gu(b, y) * dg(x, y)(c);
        });
      });
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr smat<vec<T, D, DN>, D, UP, UP>
calc_dAu(const smat<T, D, UP, UP> &gu,
         const smat<vec<T, D, DN>, D, UP, UP> &dgu, const smat<T, D, DN, DN> &A,
         const smat<vec<T, D, DN>, D, DN, DN> &dA) {
  // A^ab,c = (g^ax g^by A_xy),c
  //        = g^ax,c g^by A_xy + g^ax g^by,c A_xy + g^ax g^by A_xy,c
  return smat<vec<T, D, DN>, D, UP, UP>([&](int a, int b) ARITH_INLINE {
    return vec<T, D, DN>([&](int c) ARITH_INLINE {
      // return sum2sym([&](int x, int y) ARITH_INLINE {
      //   return dgu(a, x)(c) * gu(b, y) * A(x, y)   //
      //          + gu(a, x) * dgu(b, y)(c) * A(x, y) //
      //          + gu(a, x) * gu(b, y) * dA(x, y)(c);
      // });
      return sum<D>([&](int x) ARITH_INLINE {
        return gu(b, x) * sum<D>([&](int y) ARITH_INLINE {
                 return dgu(a, y)(c) * A(x, y)   //
                        + dgu(b, y)(c) * A(x, y) //
                        + gu(b, y) * dA(x, y)(c);
               });
      });
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr vec<smat<T, D, DN, DN>, D, DN>
calc_gammal(const smat<vec<T, D, DN>, D, DN, DN> &dg) {
  // Gammal_abc
  return vec<smat<T, D, DN, DN>, D, DN>([&](int a) ARITH_INLINE {
    return smat<T, D, DN, DN>([&](int b, int c) ARITH_INLINE {
      return (dg(a, b)(c) + dg(a, c)(b) - dg(b, c)(a)) / 2;
    });
  });
}

template <typename T, int D>
ARITH_INLINE ARITH_DEVICE ARITH_HOST constexpr vec<smat<T, D, DN, DN>, D, UP>
calc_gamma(const smat<T, D, UP, UP> &gu,
           const vec<smat<T, D, DN, DN>, D, DN> &Gammal) {
  // Gamma^a_bc
  return vec<smat<T, D, DN, DN>, D, UP>([&](int a) ARITH_INLINE {
    return smat<T, D, DN, DN>([&](int b, int c) ARITH_INLINE {
      return sum<D>([&](int x)
                        ARITH_INLINE { return gu(a, x) * Gammal(x)(b, c); });
    });
  });
}

template <typename T, int D>
ARITH_INLINE
    ARITH_DEVICE ARITH_HOST constexpr vec<vec<vec<T, D, UP>, D, DN>, D, DN>
    calc_gammalu(const smat<T, D, UP, UP> &gu,
                 const vec<smat<T, D, DN, DN>, D, DN> &Gammal) {
  // Gamma_ab^c
  return vec<vec<vec<T, D, UP>, D, DN>, D, DN>([&](int a) ARITH_INLINE {
    return vec<vec<T, D, UP>, D, DN>([&](int b) ARITH_INLINE {
      return vec<T, D, UP>([&](int c) ARITH_INLINE {
        return sum<D>([&](int x)
                          ARITH_INLINE { return Gammal(a)(b, x) * gu(x, c); });
      });
    });
  });
}

} // namespace Z4c

#endif // #ifndef PHYSICS_HXX
