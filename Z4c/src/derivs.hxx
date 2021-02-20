#ifndef DERIVS_HXX
#define DERIVS_HXX

#include "tensor.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <sstream>
#include <type_traits>

namespace Z4c {
using namespace Loop;
using namespace std;

constexpr int deriv_order = 4;

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline Z4C_INLINE T deriv1d(const T *restrict const var, const ptrdiff_t di,
                            const T dx) {
  switch (deriv_order) {
  case 2:
    return -1 / T(2) * (var[-di] - var[+di]) / dx;
  case 4:
    return (1 / T(12) * (var[-2 * di] - var[+2 * di]) //
            - 2 / T(3) * (var[-di] - var[+di])) /
           dx;
  }
  assert(0);
  return T{};
}

template <typename T>
inline Z4C_INLINE T deriv1d_upwind(const T *restrict const var,
                                   const ptrdiff_t di, const bool sign,
                                   const T dx) {
  // arXiv:1111.2177 [gr-qc], (71)
  switch (deriv_order) {
  case 2:
    if (sign)
      // +     [ 0   -1   +1    0    0]
      // + 1/2 [+1   -2   +1    0    0]

      return (1 / T(2) * var[-2 * di] //
              - 2 * var[-di]          //
              + 3 / T(2) * var[0]) /
             dx;
    else
      // +     [ 0    0   -1   +1    0  ]
      // - 1/2 [ 0    0   +1   -2   +1  ]
      //       [ 0    0   -3/2 +2   -1/2]
      return (-3 / T(2) * var[0] //
              + 2 * var[+di]     //
              - 1 / T(2) * var[+2 * di]) /
             dx;
  case 4:
    // A fourth order stencil for a first derivative, shifted by one grid point
    if (sign)
      return (-1 / T(12) * var[-3 * di] //
              + 1 / T(2) * var[-2 * di] //
              - 3 / T(2) * var[-di]     //
              + 5 / T(6) * var[0]       //
              + 1 / T(4) * var[+di]) /
             dx;
    else
      return (-1 / T(4) * var[-di]      //
              - 5 / T(6) * var[0]       //
              + 3 / T(2) * var[+di]     //
              - 1 / T(2) * var[+2 * di] //
              + 1 / T(12) * var[+3 * di]) /
             dx;
  }
  assert(0);
  return T{};
}

template <typename T>
inline Z4C_INLINE T deriv2_1d(const T *restrict const var, const ptrdiff_t di,
                              const T dx) {
  switch (deriv_order) {
  case 2:
    return ((var[-di] + var[+di]) //
            - 2 * var[0]) /
           pow2(dx);
  case 4:
    return (-1 / T(12) * (var[-2 * di] + var[+2 * di]) //
            + 4 / T(3) * (var[-di] + var[+di])         //
            - 5 / T(2) * var[0]) /
           pow2(dx);
  }
  assert(0);
  return T{};
}

template <typename T>
inline Z4C_INLINE T deriv2_2d(const T *restrict const var, const ptrdiff_t di,
                              const ptrdiff_t dj, const T dx, const T dy) {
  array<T, deriv_order + 1> arrx;
  T *const varx = &arrx[arrx.size() / 2];
  for (int j = -deriv_order / 2; j <= deriv_order / 2; ++j)
    varx[j] = deriv1d(&var[j * dj], di, dx);
  return deriv1d(varx, 1, dy);
}

template <typename T>
inline Z4C_INLINE T deriv1d_diss(const T *restrict const var,
                                 const ptrdiff_t di, const T dx) {
  switch (deriv_order) {
  case 2:
    return ((var[-2 * di] + var[+2 * di]) //
            - 4 * (var[-di] + var[+di])   //
            + 6 * var[0]) /
           dx;
  case 4:
    return ((var[-3 * di] + var[+3 * di])       //
            - 6 * (var[-2 * di] + var[+2 * di]) //
            + 15 * (var[-di] + var[+di])        //
            - 20 * var[0]) /
           dx;
  }
  assert(0);
  return T{};
}

////////////////////////////////////////////////////////////////////////////////

template <int dir, typename T>
inline Z4C_INLINE T deriv(const GF3D2<const T> &gf_, const vect<int, dim> &I,
                          const vec3<T, UP> &dx) {
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.offset(DI(dir));
  return deriv1d(&gf_(I), di, dx(dir));
}

template <int dir, typename T>
inline Z4C_INLINE T deriv_upwind(const GF3D2<const T> &gf_,
                                 const vect<int, dim> &I, const bool sign,
                                 const vec3<T, UP> &dx) {
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.offset(DI(dir));
  return deriv1d_upwind(&gf_(I), di, sign, dx(dir));
}

template <int dir1, int dir2, typename T>
inline Z4C_INLINE enable_if_t<(dir1 == dir2), T>
deriv2(const GF3D2<const T> &gf_, const vect<int, dim> &I,
       const vec3<T, UP> &dx) {
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.offset(DI(dir1));
  return deriv2_1d(&gf_(I), di, dx(dir1));
}

template <int dir1, int dir2, typename T>
inline Z4C_INLINE enable_if_t<(dir1 != dir2), T>
deriv2(const GF3D2<const T> &gf_, const vect<int, dim> &I,
       const vec3<T, UP> &dx) {
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.offset(DI(dir1));
  const ptrdiff_t dj = gf_.offset(DI(dir2));
  return deriv2_2d(&gf_(I), di, dj, dx(dir1), dx(dir2));
}

template <int dir, typename T>
inline Z4C_INLINE T deriv_diss(const GF3D2<const T> &gf_,
                               const vect<int, dim> &I, const vec3<T, UP> &dx) {
  const auto &DI = vect<int, dim>::unit;
  const ptrdiff_t di = gf_.offset(DI(dir));
  return deriv1d_diss(&gf_(I), di, dx(dir));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline Z4C_INLINE vec3<T, DN> deriv(const GF3D2<const T> &gf_,
                                    const vect<int, dim> &I,
                                    const vec3<T, UP> &dx) {
  return {
      deriv<0>(gf_, I, dx),
      deriv<1>(gf_, I, dx),
      deriv<2>(gf_, I, dx),
  };
}

template <typename T>
inline Z4C_INLINE vec3<T, DN>
deriv_upwind(const GF3D2<const T> &gf_, const vect<int, dim> &I,
             const vec3<T, UP> &dir, const vec3<T, UP> &dx) {
  return {
      deriv_upwind<0>(gf_, I, signbit(dir(0)), dx),
      deriv_upwind<1>(gf_, I, signbit(dir(1)), dx),
      deriv_upwind<2>(gf_, I, signbit(dir(2)), dx),
  };
}

template <typename T>
inline Z4C_INLINE mat3<T, DN, DN> deriv2(const GF3D2<const T> &gf_,
                                         const vect<int, dim> &I,
                                         const vec3<T, UP> &dx) {
  return {
      deriv2<0, 0>(gf_, I, dx), deriv2<0, 1>(gf_, I, dx),
      deriv2<0, 2>(gf_, I, dx), deriv2<1, 1>(gf_, I, dx),
      deriv2<1, 2>(gf_, I, dx), deriv2<2, 2>(gf_, I, dx),
  };
}

template <typename T>
inline Z4C_INLINE T diss(const GF3D2<const T> &gf_, const vect<int, dim> &I,
                         const vec3<T, UP> &dx) {
  // arXiv:gr-qc/0610128, (63), with r=2
  constexpr int diss_order = deriv_order + 2;
  constexpr int sign = diss_order % 4 == 0 ? -1 : +1;
  return sign / T(pown(2, deriv_order + 2)) *
         (deriv_diss<0>(gf_, I, dx)   //
          + deriv_diss<1>(gf_, I, dx) //
          + deriv_diss<2>(gf_, I, dx));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH, const GF3D2<const T> &gf0_,
            const GF3D2<T> &gf_, const vec3<GF3D2<T>, DN> &dgf_) {
  DECLARE_CCTK_ARGUMENTS;

  const vec3<CCTK_REAL, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    const auto val = gf0_(p.I);
    gf_(p.I) = val;
    const auto dval = deriv(gf0_, p.I, dx);
    for (int a = 0; a < 3; ++a)
      dgf_(a)(p.I) = dval(a);
  });
}

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH, const GF3D2<const T> &gf0_,
             const GF3D2<T> &gf_, const vec3<GF3D2<T>, DN> &dgf_,
             const mat3<GF3D2<T>, DN, DN> &ddgf_) {
  DECLARE_CCTK_ARGUMENTS;

  const vec3<CCTK_REAL, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    const auto val = gf0_(p.I);
    gf_(p.I) = val;
    const auto dval = deriv(gf0_, p.I, dx);
    for (int a = 0; a < 3; ++a)
      dgf_(a)(p.I) = dval(a);
    const auto ddval = deriv2(gf0_, p.I, dx);
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        ddgf_(a, b)(p.I) = ddval(a, b);
  });
}

template <typename T, dnup_t dnup>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH,
            const vec3<GF3D2<const T>, dnup> &gf0_,
            const vec3<GF3D2<T>, dnup> &gf_,
            const vec3<vec3<GF3D2<T>, DN>, dnup> &dgf_) {
  for (int a = 0; a < 3; ++a)
    calc_derivs(cctkGH, gf0_(a), gf_(a), dgf_(a));
}

template <typename T, dnup_t dnup>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const cGH *restrict const cctkGH, const vec3<GF3D2<const T>, dnup> &gf0_,
    const vec3<GF3D2<T>, dnup> &gf_, const vec3<vec3<GF3D2<T>, DN>, dnup> &dgf_,
    const vec3<mat3<GF3D2<T>, DN, DN>, dnup> &ddgf_) {
  for (int a = 0; a < 3; ++a)
    calc_derivs2(cctkGH, gf0_(a), gf_(a), dgf_(a), ddgf_(a));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH,
            const mat3<GF3D2<const T>, dnup1, dnup2> &gf0_,
            const mat3<GF3D2<T>, dnup1, dnup2> &gf_,
            const mat3<vec3<GF3D2<T>, DN>, dnup1, dnup2> &dgf_) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs(cctkGH, gf0_(a, b), gf_(a, b), dgf_(a, b));
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH,
             const mat3<GF3D2<const T>, dnup1, dnup2> &gf0_,
             const mat3<GF3D2<T>, dnup1, dnup2> &gf_,
             const mat3<vec3<GF3D2<T>, DN>, dnup1, dnup2> &dgf_,
             const mat3<mat3<GF3D2<T>, DN, DN>, dnup1, dnup2> &ddgf_) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs2(cctkGH, gf0_(a, b), gf_(a, b), dgf_(a, b), ddgf_(a, b));
}

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
apply_upwind_diss(const cGH *restrict const cctkGH, const GF3D2<const T> &gf_,
                  const vec3<GF3D2<const T>, UP> &gf_betaG_,
                  const GF3D2<T> &gf_rhs_) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const vec3<CCTK_REAL, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  if (epsdiss == 0) {

    loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
      const vec3<CCTK_REAL, UP> betaG = gf_betaG_(p.I);
      const vec3<CCTK_REAL, DN> dgf_upwind(deriv_upwind(gf_, p.I, betaG, dx));
      gf_rhs_(p.I) += sum1([&](int x) { return betaG(x) * dgf_upwind(x); });
    });

  } else {

    loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
      const vec3<CCTK_REAL, UP> betaG = gf_betaG_(p.I);
      const vec3<CCTK_REAL, DN> dgf_upwind(deriv_upwind(gf_, p.I, betaG, dx));
      gf_rhs_(p.I) += sum1([&](int x) { return betaG(x) * dgf_upwind(x); }) //
                      + epsdiss * diss(gf_, p.I, dx);
    });
  }
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH, const GF3D2<const T> &gf1,
            const GF3D5<T> &gf0, const vec3<GF3D5<T>, DN> &dgf0,
            const GF3D5layout &layout0) {
  DECLARE_CCTK_ARGUMENTS;

  const vec3<CCTK_REAL, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    const auto val = gf1(p.I);
    gf0(layout0, p.I) = val;
    const auto dval = deriv(gf1, p.I, dx);
    for (int a = 0; a < 3; ++a)
      dgf0(a)(layout0, p.I) = dval(a);
  });
}

template <typename T>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH, const GF3D2<const T> &gf1,
             const GF3D5<T> &gf0, const vec3<GF3D5<T>, DN> &dgf0,
             const mat3<GF3D5<T>, DN, DN> &ddgf0, const GF3D5layout &layout0) {
  DECLARE_CCTK_ARGUMENTS;

  const vec3<CCTK_REAL, UP> dx([&](int a) { return CCTK_DELTA_SPACE(a); });

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) Z4C_INLINE {
    const auto val = gf1(p.I);
    gf0(layout0, p.I) = val;
    const auto dval = deriv(gf1, p.I, dx);
    for (int a = 0; a < 3; ++a)
      dgf0(a)(layout0, p.I) = dval(a);
    const auto ddval = deriv2(gf1, p.I, dx);
    for (int a = 0; a < 3; ++a)
      for (int b = a; b < 3; ++b)
        ddgf0(a, b)(layout0, p.I) = ddval(a, b);
  });
}

template <typename T, dnup_t dnup>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs(
    const cGH *restrict const cctkGH, const vec3<GF3D2<const T>, dnup> &gf0_,
    const vec3<GF3D5<T>, dnup> &gf_, const vec3<vec3<GF3D5<T>, DN>, dnup> &dgf_,
    const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    calc_derivs(cctkGH, gf0_(a), gf_(a), dgf_(a), layout);
}

template <typename T, dnup_t dnup>
CCTK_ATTRIBUTE_NOINLINE void calc_derivs2(
    const cGH *restrict const cctkGH, const vec3<GF3D2<const T>, dnup> &gf0_,
    const vec3<GF3D5<T>, dnup> &gf_, const vec3<vec3<GF3D5<T>, DN>, dnup> &dgf_,
    const vec3<mat3<GF3D5<T>, DN, DN>, dnup> &ddgf_,
    const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    calc_derivs2(cctkGH, gf0_(a), gf_(a), dgf_(a), ddgf_(a), layout);
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs(const cGH *restrict const cctkGH,
            const mat3<GF3D2<const T>, dnup1, dnup2> &gf0_,
            const mat3<GF3D5<T>, dnup1, dnup2> &gf_,
            const mat3<vec3<GF3D5<T>, DN>, dnup1, dnup2> &dgf_,
            const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs(cctkGH, gf0_(a, b), gf_(a, b), dgf_(a, b), layout);
}

template <typename T, dnup_t dnup1, dnup_t dnup2>
CCTK_ATTRIBUTE_NOINLINE void
calc_derivs2(const cGH *restrict const cctkGH,
             const mat3<GF3D2<const T>, dnup1, dnup2> &gf0_,
             const mat3<GF3D5<T>, dnup1, dnup2> &gf_,
             const mat3<vec3<GF3D5<T>, DN>, dnup1, dnup2> &dgf_,
             const mat3<mat3<GF3D5<T>, DN, DN>, dnup1, dnup2> &ddgf_,
             const GF3D5layout &layout) {
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b)
      calc_derivs2(cctkGH, gf0_(a, b), gf_(a, b), dgf_(a, b), ddgf_(a, b),
                   layout);
}

} // namespace Z4c

#endif // #ifndef DERIVS_HXX
