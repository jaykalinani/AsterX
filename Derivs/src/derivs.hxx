#ifndef DERIVS_HXX
#define DERIVS_HXX

#include <defs.hxx>
#include <div.hxx>
#include <loop_device.hxx>
#include <mat.hxx>
#include <simd.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <type_traits>

namespace Derivs {

////////////////////////////////////////////////////////////////////////////////

namespace detail {
using namespace Arith;
using namespace Loop;

// Pointwise one-dimensional derivative operators

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST simd<T>
copy1d(const simdl<T> &mask, const T *restrict const var) {
  return maskz_loadu(mask, var);
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, R>
    deriv1d(const TS var, const T dx) {
  constexpr T c1 = 1 / T(2);
  return c1 * (var(1) - var(-1)) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, R>
    deriv1d(const TS var, const T dx) {
  constexpr T c1 = 2 / T(3);
  constexpr T c2 = -1 / T(12);
  return (c2 * (var(2) - var(-2)) + c1 * (var(1) - var(-1))) / dx;
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, simd<T> >
    deriv1d_upwind(const simdl<T> &mask, const T *restrict const var,
                   const std::ptrdiff_t di, const simd<T> &vel, const T dx) {
  // arXiv:1111.2177 [gr-qc], (71)

  // if (sign)
  //   // +     [ 0   -1   +1    0    0]
  //   // + 1/2 [+1   -2   +1    0    0]
  //   //       [+1/2 -2   +3/2  0    0]
  //   return (1 / T(2) * var[-2 * di] //
  //           - 2 * var[-di]          //
  //           + 3 / T(2) * var[0]) /
  //          dx;
  // else
  //   // +     [ 0    0   -1   +1    0  ]
  //   // - 1/2 [ 0    0   +1   -2   +1  ]
  //   //       [ 0    0   -3/2 +2   -1/2]
  //   return (-3 / T(2) * var[0] //
  //           + 2 * var[+di]     //
  //           - 1 / T(2) * var[+2 * di]) /
  //          dx;

  // + 1/2 [+1/2 -2   +3/2  0    0  ]
  // + 1/2 [ 0    0   -3/2 +2   -1/2]
  //       [+1/4 -1    0   +1   -1/4]
  constexpr T c1s = 1;
  constexpr T c2s = -1 / T(4);
  const simd<T> symm =
      c2s * (maskz_loadu(mask, &var[2 * di]) -
             maskz_loadu(mask, &var[-2 * di])) //
      + c1s(maskz_loadu(mask, &var[di]) - maskz_loadu(mask, &var[-di]));
  // + 1/2 [+1/2 -2   +3/2  0    0  ]
  // - 1/2 [ 0    0   -3/2 +2   -1/2]
  //       [+1/4 -1   +3/2 -1   +1/4]
  constexpr T c0a = 3 / T(2);
  constexpr T c1a = -1;
  constexpr T c2a = 1 / T(4);
  const simd<T> anti =
      c2a * (maskz_loadu(mask, &var[2 * di]) +
             maskz_loadu(mask, &var[-2 * di]))                          //
      + c1a(maskz_loadu(mask, &var[di]) + maskz_loadu(mask, &var[-di])) //
      + c0a * maskz_loadu(mask, &var[0]);
  using std::fabs;
  return (vel * symm - fabs(vel) * anti) / dx;
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, simd<T> >
    deriv1d_upwind(const simdl<T> &mask, const T *restrict const var,
                   const std::ptrdiff_t di, const simd<T> &vel, const T dx) {
  // arXiv:1111.2177 [gr-qc], (71)

  // A fourth order stencil for a first derivative, shifted by one grid point

  // if (sign)
  //   return (-1 / T(12) * var[-3 * di] //
  //           + 1 / T(2) * var[-2 * di] //
  //           - 3 / T(2) * var[-di]     //
  //           + 5 / T(6) * var[0]       //
  //           + 1 / T(4) * var[+di]) /
  //          dx;
  // else
  //   return (-1 / T(4) * var[-di]      //
  //           - 5 / T(6) * var[0]       //
  //           + 3 / T(2) * var[+di]     //
  //           - 1 / T(2) * var[+2 * di] //
  //           + 1 / T(12) * var[+3 * di]) /
  //          dx;

  constexpr T c1s = 7 / T(8);
  constexpr T c2s = -1 / T(4);
  constexpr T c3s = 1 / T(24);
  const simd<T> symm =
      c3s * (maskz_loadu(mask, &var[3 * di]) -
             maskz_loadu(mask, &var[-3 * di])) //
      + c2s * (maskz_loadu(mask, &var[2 * di]) -
               maskz_loadu(mask, &var[-2 * di])) //
      + c1s * (maskz_loadu(mask, &var[di]) - maskz_loadu(mask, &var[-di]));
  constexpr T c0a = 5 / T(6);
  constexpr T c1a = -5 / T(8);
  constexpr T c2a = 1 / T(4);
  constexpr T c3a = -1 / T(24);
  const simd<T> anti =
      c3a * (maskz_loadu(mask, &var[3 * di]) +
             maskz_loadu(mask, &var[-3 * di])) //
      + c2a * (maskz_loadu(mask, &var[2 * di]) +
               maskz_loadu(mask, &var[-2 * di]))                           //
      + c1a * (maskz_loadu(mask, &var[di]) + maskz_loadu(mask, &var[-di])) //
      + c0a * maskz_loadu(mask, &var[0]);
  using std::fabs;
  return (vel * symm - fabs(vel) * anti) / dx;
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, R>
    deriv2_1d(const TS var, const T dx) {
  constexpr T c0 = -2;
  constexpr T c1 = 1;
  return (c1 * (var(-1) + var(1)) + c0 * var(0)) / pow2(dx);
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, R>
    deriv2_1d(const TS var, const T dx) {
  constexpr T c0 = -5 / T(2);
  constexpr T c1 = 4 / T(3);
  constexpr T c2 = -1 / T(12);
  return (c2 * (var(-2) + var(2)) + c1 * (var(-1) + var(1)) + c0 * var(0)) /
         pow2(dx);
}

template <int deriv_order, typename T, typename TS,
          typename R = std::result_of_t<TS(int, int)> >
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
deriv2_2d(const TS var, const T dx, const T dy) {
  // We assume that the x-direction might be special since it might
  // be SIMD-vectorized. We assume that the y-direction is not
  // SIMD-vectorized.

  // Calculate y-derivative first.
  // (If we wanted to calculate the x-derivative first, then we would
  // be difficult to determine the extent of the `n`-loop, since it
  // needs to include enough room for the SIMD y-derivative later.)
  static_assert(sizeof(R) % sizeof(T) == 0, "");
  constexpr std::ptrdiff_t vsize = sizeof(R) / sizeof(T);
  constexpr std::ptrdiff_t ndyvars =
      align_ceil(std::ptrdiff_t(2 * deriv_order + 1), vsize);
  std::array<R, ndyvars> dyvar;
#ifdef CCTK_DEBUG
  for (std::ptrdiff_t n = 0; n < ndyvars; ++n)
    dyvar[n] = Arith::nan<T>()();
#endif
  for (std::ptrdiff_t n = 0; n < ndyvars; ++n) {
    std::ptrdiff_t di = vsize * n - deriv_order;
    if (vsize == 1 && di == 0)
      continue;
    dyvar[n] = deriv1d<deriv_order>(
        [&](int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE { return var(di, dj); }, dy);
  }

  // Calculate x-derivative next
  return deriv1d<deriv_order>(
      [&](int di)
          CCTK_ATTRIBUTE_ALWAYS_INLINE { return dyvar[di + deriv_order]; },
      dx);
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 2, simd<T> >
    diss1d(const simdl<T> &mask, const T *restrict const var,
           const std::ptrdiff_t di, const T dx) {
  constexpr T c0 = 6;
  constexpr T c1 = -4;
  constexpr T c2 = 1;
  return (c2 * (maskz_loadu(mask, &var[2 * di]) +
                maskz_loadu(mask, &var[-2 * di]))                             //
          + c1 * (maskz_loadu(mask, &var[di]) + maskz_loadu(mask, &var[-di])) //
          + c0 * maskz_loadu(mask, &var[0])) /
         dx;
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST std::enable_if_t<deriv_order == 4, simd<T> >
    diss1d(const simdl<T> &mask, const T *restrict const var,
           const std::ptrdiff_t di, const T dx) {
  constexpr T c0 = -20;
  constexpr T c1 = 15;
  constexpr T c2 = -6;
  constexpr T c3 = 1;
  return (c3 * (maskz_loadu(mask, &var[3 * di]) +
                maskz_loadu(mask, &var[-3 * di])) //
          + c2 * *(maskz_loadu(mask, &var[2 * di]) +
                   maskz_loadu(mask, &var[-2 * di]))                          //
          + c1 * (maskz_loadu(mask, &var[di]) + maskz_loadu(mask, &var[-di])) //
          + c0 * maskz_loadu(mask, &var[0])) /
         dx;
}

} // namespace detail

////////////////////////////////////////////////////////////////////////////////

// Pointwise multi-dimensional derivative operators

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST Arith::vec<Arith::simd<T>, Loop::dim, Arith::DN>
    calc_deriv(const Loop::GF3D5<const T> &gf, const Arith::simdl<T> &mask,
               const Loop::GF3D5layout &layout,
               const Arith::vect<int, Loop::dim> &I,
               const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(layout, I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      layout.delta(1, 0, 0),
      layout.delta(0, 1, 0),
      layout.delta(0, 0, 1),
  };
  return {
      detail::deriv1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[0]]);
          },
          dx[0]),
      detail::deriv1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[1]]);
          },
          dx[1]),
      detail::deriv1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[2]]);
          },
          dx[2]),
  };
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST Arith::vec<T, Loop::dim, Arith::DN>
    calc_deriv(const Loop::GF3D5<const T> &gf, const Loop::GF3D5layout &layout,
               const Arith::vect<int, Loop::dim> &I,
               const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(layout, I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      layout.delta(1, 0, 0),
      layout.delta(0, 1, 0),
      layout.delta(0, 0, 1),
  };
  return {
      detail::deriv1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[0]]; },
          dx[0]),
      detail::deriv1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[1]]; },
          dx[1]),
      detail::deriv1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[2]]; },
          dx[2]),
  };
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
    Arith::smat<Arith::simd<T>, Loop::dim, Arith::DN, Arith::DN>
    calc_deriv2(const Loop::GF3D5<const T> &gf, const Arith::simdl<T> &mask,
                const Loop::GF3D5layout &layout,
                const Arith::vect<int, Loop::dim> &I,
                const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(layout, I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      layout.delta(1, 0, 0),
      layout.delta(0, 1, 0),
      layout.delta(0, 0, 1),
  };
  return {
      detail::deriv2_1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[0]]);
          },
          dx[0]),
      detail::deriv2_2d<deriv_order>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[0] + dj * offsets[1]]);
          },
          dx[0], dx[1]),
      detail::deriv2_2d<deriv_order>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[0] + dj * offsets[2]]);
          },
          dx[0], dx[2]),
      detail::deriv2_1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[1]]);
          },
          dx[1]),
      detail::deriv2_2d<deriv_order>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[1] + dj * offsets[2]]);
          },
          dx[1], dx[2]),
      detail::deriv2_1d<deriv_order>(
          [&](int di) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return maskz_loadu(mask, &ptr[di * offsets[2]]);
          },
          dx[2]),
  };
}

template <int deriv_order, typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE
    CCTK_DEVICE CCTK_HOST Arith::smat<T, Loop::dim, Arith::DN, Arith::DN>
    calc_deriv2(const Loop::GF3D5<const T> &gf, const Loop::GF3D5layout &layout,
                const Arith::vect<int, Loop::dim> &I,
                const Arith::vect<T, Loop::dim> &dx) {
  using namespace Arith;
  using namespace Loop;
  // We use explicit index calculations to avoid unnecessary integer
  // multiplications
  const T *restrict const ptr = &gf(layout, I);
  const std::array<std::ptrdiff_t, Loop::dim> offsets{
      layout.delta(1, 0, 0),
      layout.delta(0, 1, 0),
      layout.delta(0, 0, 1),
  };
  return {
      detail::deriv2_1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[0]]; },
          dx[0]),
      detail::deriv2_2d<deriv_order>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return ptr[di * offsets[0] + dj * offsets[1]];
          },
          dx[0], dx[1]),
      detail::deriv2_2d<deriv_order>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return ptr[di * offsets[0] + dj * offsets[2]];
          },
          dx[0], dx[2]),
      detail::deriv2_1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[1]]; },
          dx[1]),
      detail::deriv2_2d<deriv_order>(
          [&](int di, int dj) CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return ptr[di * offsets[1] + dj * offsets[2]];
          },
          dx[1], dx[2]),
      detail::deriv2_1d<deriv_order>(
          [&](int di)
              CCTK_ATTRIBUTE_ALWAYS_INLINE { return ptr[di * offsets[2]]; },
          dx[2]),
  };
}

////////////////////////////////////////////////////////////////////////////////

// Tile-based  multi-dimensional operators

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const Arith::vec<Loop::GF3D5<T>, Loop::dim, Arith::DN> &dgf,
            const Loop::GridDescBaseDevice &grid,
            const Loop::GF3D5<const T> &gf, const Loop::GF3D5layout layout,
            const Arith::vect<T, Loop::dim> dx, const int deriv_order);

template <int CI, int CJ, int CK, typename T>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const Arith::vec<Loop::GF3D5<T>, Loop::dim, Arith::DN> &dgf,
    const Arith::smat<Loop::GF3D5<T>, Loop::dim, Arith::DN, Arith::DN> &ddgf,
    const Loop::GridDescBaseDevice &grid, const Loop::GF3D5<const T> &gf,
    const Loop::GF3D5layout layout, const Arith::vect<T, Loop::dim> dx,
    const int deriv_order);

} // namespace Derivs

#endif // #ifndef DERIVS_HXX
