#ifndef SIMD_HXX
#define SIMD_HXX

#ifdef __CUDACC__
#define CCTK_CUDA __device__ __host__
#else
#define CCTK_CUDA
#endif

#include <fixmath.hxx>
#include <vect.hxx>

#include <nsimd/nsimd-all.hpp>
#undef vec // This should arguably not be defined in C++

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <utility>
#include <type_traits>

namespace Arith {

template <typename T> struct simd;
template <typename T> struct simdl;

namespace detail {
template <typename T> struct f2i;
template <typename T> using f2i_t = typename f2i<T>::type;
template <> struct f2i<f32> { typedef i32 type; };
template <> struct f2i<f64> { typedef i64 type; };

template <typename T> struct f2u;
template <typename T> using f2u_t = typename f2u<T>::type;
template <> struct f2u<f32> { typedef u32 type; };
template <> struct f2u<f64> { typedef u64 type; };
} // namespace detail

template <typename T> struct simd {
  typedef T value_type;
#ifndef SIMD_CPU
  nsimd::pack<T> elts;
#else
  T elts;
#endif

  simd(const simd &) = default;
  simd(simd &&) = default;
  simd &operator=(const simd &) = default;
  simd &operator=(simd &&) = default;

  CCTK_CUDA simd() {}
  CCTK_CUDA simd(const T &a) : elts(a) {}
#ifndef SIMD_CPU
  CCTK_CUDA simd(const nsimd::pack<T> &elts) : elts(elts) {}
#endif

  CCTK_CUDA constexpr std::size_t size() const {
#ifndef SIMD_CPU
    return sizeof(nsimd::pack<T>) / sizeof(T);
#else
    return 1;
#endif
  }

  CCTK_CUDA simdl<T> operator!() const { return !elts; }
  CCTK_CUDA simd operator~() const { return ~elts; }

  CCTK_CUDA simd operator+() const { return +elts; }
  CCTK_CUDA simd operator-() const { return -elts; }

  friend CCTK_CUDA simd operator&(const simd &x, const simd &y) {
    return x.elts & y.elts;
  }
  friend CCTK_CUDA simd operator|(const simd &x, const simd &y) {
    return x.elts | y.elts;
  }
  friend CCTK_CUDA simd operator^(const simd &x, const simd &y) {
    return x.elts ^ y.elts;
  }
  friend CCTK_CUDA simd operator+(const simd &x, const simd &y) {
    return x.elts + y.elts;
  }
  friend CCTK_CUDA simd operator-(const simd &x, const simd &y) {
    return x.elts - y.elts;
  }
  friend CCTK_CUDA simd operator*(const simd &x, const simd &y) {
    return x.elts * y.elts;
  }
  friend CCTK_CUDA simd operator/(const simd &x, const simd &y) {
    return x.elts / y.elts;
  }
  friend CCTK_CUDA simd operator%(const simd &x, const simd &y) {
    return x.elts % y.elts;
  }

  friend CCTK_CUDA simd operator&(const T &a, const simd &y) {
    return a & y.elts;
  }
  friend CCTK_CUDA simd operator|(const T &a, const simd &y) {
    return a | y.elts;
  }
  friend CCTK_CUDA simd operator^(const T &a, const simd &y) {
    return a ^ y.elts;
  }
  friend CCTK_CUDA simd operator+(const T &a, const simd &y) {
    return a + y.elts;
  }
  friend CCTK_CUDA simd operator-(const T &a, const simd &y) {
    return a - y.elts;
  }
  friend CCTK_CUDA simd operator*(const T &a, const simd &y) {
    return a * y.elts;
  }
  friend CCTK_CUDA simd operator/(const T &a, const simd &y) {
    return a / y.elts;
  }
  friend CCTK_CUDA simd operator%(const T &a, const simd &y) {
    return a % y.elts;
  }

  friend CCTK_CUDA simd operator&(const simd &x, const T &b) {
    return x.elts & b;
  }
  friend CCTK_CUDA simd operator|(const simd &x, const T &b) {
    return x.elts | b;
  }
  friend CCTK_CUDA simd operator^(const simd &x, const T &b) {
    return x.elts ^ b;
  }
  friend CCTK_CUDA simd operator+(const simd &x, const T &b) {
    return x.elts + b;
  }
  friend CCTK_CUDA simd operator-(const simd &x, const T &b) {
    return x.elts - b;
  }
  friend CCTK_CUDA simd operator*(const simd &x, const T &b) {
    return x.elts * b;
  }
  friend CCTK_CUDA simd operator/(const simd &x, const T &b) {
    return x.elts / b;
  }
  friend CCTK_CUDA simd operator%(const simd &x, const T &b) {
    return x.elts % b;
  }

  CCTK_CUDA simd &operator&=(const simd &x) { return *this = *this & x; }
  CCTK_CUDA simd &operator|=(const simd &x) { return *this = *this | x; }
  CCTK_CUDA simd &operator^=(const simd &x) { return *this = *this ^ x; }
  CCTK_CUDA simd &operator+=(const simd &x) { return *this = *this + x; }
  CCTK_CUDA simd &operator-=(const simd &x) { return *this = *this - x; }
  CCTK_CUDA simd &operator*=(const simd &x) { return *this = *this * x; }
  CCTK_CUDA simd &operator/=(const simd &x) { return *this = *this / x; }
  CCTK_CUDA simd &operator%=(const simd &x) { return *this = *this % x; }

  CCTK_CUDA simd &operator&=(const T &a) { return *this = *this & a; }
  CCTK_CUDA simd &operator|=(const T &a) { return *this = *this | a; }
  CCTK_CUDA simd &operator^=(const T &a) { return *this = *this ^ a; }
  CCTK_CUDA simd &operator+=(const T &a) { return *this = *this + a; }
  CCTK_CUDA simd &operator-=(const T &a) { return *this = *this - a; }
  CCTK_CUDA simd &operator*=(const T &a) { return *this = *this * a; }
  CCTK_CUDA simd &operator/=(const T &a) { return *this = *this / a; }
  CCTK_CUDA simd &operator%=(const T &a) { return *this = *this % a; }

  friend CCTK_CUDA simdl<T> operator==(const simd &x, const simd &y) {
    return x.elts == y.elts;
  }
  friend CCTK_CUDA simdl<T> operator!=(const simd &x, const simd &y) {
    return x.elts != y.elts;
  }
  friend CCTK_CUDA simdl<T> operator<(const simd &x, const simd &y) {
    return x.elts < y.elts;
  }
  friend CCTK_CUDA simdl<T> operator>(const simd &x, const simd &y) {
    return x.elts > y.elts;
  }
  friend CCTK_CUDA simdl<T> operator<=(const simd &x, const simd &y) {
    return x.elts <= y.elts;
  }
  friend CCTK_CUDA simdl<T> operator>=(const simd &x, const simd &y) {
    return x.elts >= y.elts;
  }

  friend CCTK_CUDA simdl<T> operator==(const T &a, const simd &y) {
    return a == y.elts;
  }
  friend CCTK_CUDA simdl<T> operator!=(const T &a, const simd &y) {
    return a != y.elts;
  }
  friend CCTK_CUDA simdl<T> operator<(const T &a, const simd &y) {
    return a < y.elts;
  }
  friend CCTK_CUDA simdl<T> operator>(const T &a, const simd &y) {
    return a > y.elts;
  }
  friend CCTK_CUDA simdl<T> operator<=(const T &a, const simd &y) {
    return a <= y.elts;
  }
  friend CCTK_CUDA simdl<T> operator>=(const T &a, const simd &y) {
    return a >= y.elts;
  }

  friend CCTK_CUDA simdl<T> operator==(const simd &x, const T &b) {
    return x.elts == b;
  }
  friend CCTK_CUDA simdl<T> operator!=(const simd &x, const T &b) {
    return x.elts != b;
  }
  friend CCTK_CUDA simdl<T> operator<(const simd &x, const T &b) {
    return x.elts < b;
  }
  friend CCTK_CUDA simdl<T> operator>(const simd &x, const T &b) {
    return x.elts > b;
  }
  friend CCTK_CUDA simdl<T> operator<=(const simd &x, const T &b) {
    return x.elts <= b;
  }
  friend CCTK_CUDA simdl<T> operator>=(const simd &x, const T &b) {
    return x.elts >= b;
  }

  friend CCTK_CUDA simd abs(const simd &x) {
    using std::abs;
    return abs(x.elts);
  }
  friend CCTK_CUDA simd fabs(const simd &x) {
    using std::abs;
    return abs(x.elts);
  }
  friend CCTK_CUDA simdl<T> signbit(const simd &x) {
#ifndef SIMD_CPU
    typedef detail::f2u_t<T> U;
    constexpr T signmask = scalar_reinterpret(T{}, U(1) << (8 * sizeof(U) - 1));
    return to_logical(x & signmask);
#else
    using std::signbit;
    return signbit(x.elts);
#endif
  }
  friend CCTK_CUDA simd sqrt(const simd &x) {
    using std::sqrt;
    return sqrt(x.elts);
  }
  friend CCTK_CUDA simdl<T> to_logical(const simd &x) {
    return to_logical(x.elts);
  }

  friend CCTK_CUDA simd copysign(const simd &x, const simd &y) {
#ifndef SIMD_CPU
    typedef detail::f2u_t<T> U;
    constexpr T signmask = scalar_reinterpret(T{}, U(1) << (8 * sizeof(U) - 1));
    return (x & ~signmask) | (y & signmask);
#else
    using std::copysign;
    return copysign(x.elts, y.elts);
#endif
  }
  friend CCTK_CUDA simd flipsign(const simd &x, const simd &y) {
#ifndef SIMD_CPU
    typedef detail::f2u_t<T> U;
    constexpr T signmask = scalar_reinterpret(T{}, U(1) << (8 * sizeof(U) - 1));
    return x ^ (y & signmask);
#else
    return copysign(1, y) * x;
#endif
  }
  friend CCTK_CUDA simd fmax(const simd &x, const simd &y) {
    using std::max;
    return max(x.elts, y.elts);
  }
  friend CCTK_CUDA simd fmin(const simd &x, const simd &y) {
    using std::min;
    return min(x.elts, y.elts);
  }
  friend CCTK_CUDA simd max(const simd &x, const simd &y) {
    using std::max;
    return max(x.elts, y.elts);
  }
  friend CCTK_CUDA simd min(const simd &x, const simd &y) {
    using std::min;
    return min(x.elts, y.elts);
  }

  friend CCTK_CUDA simd fmax(const T &a, const simd &y) {
    using std::max;
    return max(a, y.elts);
  }
  friend CCTK_CUDA simd fmin(const T &a, const simd &y) {
    using std::min;
    return min(a, y.elts);
  }
  friend CCTK_CUDA simd max(const T &a, const simd &y) {
    using std::max;
    return max(a, y.elts);
  }
  friend CCTK_CUDA simd min(const T &a, const simd &y) {
    using std::min;
    return min(a, y.elts);
  }

  friend CCTK_CUDA simd fmax(const simd &x, const T &b) {
    using std::max;
    return max(x.elts, b);
  }
  friend CCTK_CUDA simd fmin(const simd &x, const T &b) {
    using std::min;
    return min(x.elts, b);
  }
  friend CCTK_CUDA simd max(const simd &x, const T &b) {
    using std::max;
    return max(x.elts, b);
  }
  friend CCTK_CUDA simd min(const simd &x, const T &b) {
    using std::min;
    return min(x.elts, b);
  }

  friend CCTK_CUDA void storea(T *ptr, const simd &x) {
#ifndef SIMD_CPU
    nsimd::storea(ptr, x.elts);
#else
    *ptr = x.elts;
#endif
  }
  friend CCTK_CUDA void storeu(T *ptr, const simd &x) {
#ifndef SIMD_CPU
    nsimd::storeu(ptr, x.elts);
#else
    *ptr = x.elts;
#endif
  }
  friend CCTK_CUDA void mask_storea(const simdl<T> &mask, T *ptr,
                                    const simd &x) {
#ifndef SIMD_CPU
    nsimd::mask_storea(mask.elts, ptr, x.elts);
#else
    if (mask.elts)
      *ptr = x.elts;
#endif
  }
  friend CCTK_CUDA void mask_storeu(const simdl<T> &mask, T *ptr,
                                    const simd &x) {
#ifndef SIMD_CPU
    nsimd::mask_storeu(mask.elts, ptr, x.elts);
#else
    if (mask.elts)
      *ptr = x.elts;
#endif
  }

  friend ostream &operator<<(ostream &os, const simd &x) {
    return os << "simd";
  }
};

} // namespace Arith
namespace std {
template <typename T>
struct tuple_size<Arith::simd<T> >
    : std::integral_constant<std::size_t, sizeof(Arith::simd<T>) / sizeof(T)> {
};
} // namespace std
namespace Arith {

template <typename T> struct zero<Arith::simd<T> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE Arith::simd<T> operator()() const {
    return 0;
  }
};

template <typename VT, typename T = typename VT::value_type>
CCTK_CUDA inline simd<T> iota() {
#ifndef SIMD_CPU
  return nsimd::iota<nsimd::pack<T> >();
#else
  return 0;
#endif
}

template <typename VT, typename T = typename VT::value_type>
CCTK_CUDA inline simdl<T> mask_for_loop_tail(const int i, const int n) {
#ifndef SIMD_CPU
  return nsimd::mask_for_loop_tail<nsimd::packl<T> >(i, n);
#else
  return i < n;
#endif
}

template <typename VT, typename T = typename VT::value_type>
CCTK_CUDA inline simd<T> loada(const T *ptr) {
#ifndef SIMD_CPU
  return nsimd::loada<nsimd::pack<T> >(ptr);
#else
  return *ptr;
#endif
}

template <typename VT, typename T = typename VT::value_type>
CCTK_CUDA inline simd<T> loadu(const T *ptr) {
#ifndef SIMD_CPU
  return nsimd::loadu<nsimd::pack<T> >(ptr);
#else
  return *ptr;
#endif
}

template <typename T>
CCTK_CUDA inline simd<T> maskz_loada(const simdl<T> &mask, const T *ptr) {
#ifndef SIMD_CPU
  return nsimd::maskz_loada(mask.elts, ptr);
#else
  return mask.elts ? *ptr : 0;
#endif
}

template <typename T>
CCTK_CUDA inline simd<T> maskz_loadu(const simdl<T> &mask, const T *ptr) {
#ifndef SIMD_CPU
  return nsimd::maskz_loadu(mask.elts, ptr);
#else
  return mask.elts ? *ptr : 0;
#endif
}

template <typename T>
CCTK_CUDA inline simd<T> masko_loada(const simdl<T> &mask, const T *ptr,
                                     const simd<T> &other) {
#ifndef SIMD_CPU
  return nsimd::masko_loada(mask.elts, ptr, other.elts);
#else
  return mask.elts ? *ptr : other.elts;
#endif
}

template <typename T>
CCTK_CUDA inline simd<T> masko_loadu(const simdl<T> &mask, const T *ptr,
                                     const simd<T> &other) {
#ifndef SIMD_CPU
  return nsimd::masko_loadu(mask.elts, ptr, other.elts);
#else
  return mask.elts ? *ptr : other.elts;
#endif
}

template <typename T, typename U,
          enable_if_t<is_convertible_v<T, U> > * = nullptr>
CCTK_CUDA inline simd<T> masko_loada(const simdl<T> &mask, const T *ptr,
                                     const U &other) {
#ifndef SIMD_CPU
  return masko_loada(mask, ptr, simd<T>(other));
#else
  return mask.elts ? *ptr : other;
#endif
}

template <typename T, typename U,
          enable_if_t<is_convertible_v<T, U> > * = nullptr>
CCTK_CUDA inline simd<T> masko_loadu(const simdl<T> &mask, const T *ptr,
                                     const U &other) {
#ifndef SIMD_CPU
  return masko_loadu(mask, ptr, simd<T>(other));
#else
  return mask.elts ? *ptr : other;
#endif
}

template <typename T> CCTK_CUDA inline simd<T> cbrt(const simd<T> &x) {
  // alignas(alignof(simd<T>)) T xarr[x.size()];
  // storea(xarr, x);
  // alignas(alignof(simd<T>)) T yarr[x.size()];
  // for (std::size_t n = 0; n < x.size(); ++n) {
  //   using std::cbrt;
  //   yarr[n] = cbrt(xarr[n]);
  // }
  // const simd<T> y = loada<simd<T> >(yarr);
  // return y;
  constexpr size_t vsize = tuple_size_v<simd<T> >;
  T xarr[vsize];
  storeu(xarr, x);
  T yarr[vsize];
  for (std::size_t n = 0; n < vsize; ++n) {
    using std::cbrt;
    yarr[n] = cbrt(xarr[n]);
  }
  const simd<T> y = loadu<simd<T> >(yarr);
  return y;
}

template <typename T> CCTK_CUDA inline simd<T> sin(const simd<T> &x) {
  // alignas(alignof(simd<T>)) T xarr[x.size()];
  // storea(xarr, x);
  // alignas(alignof(simd<T>)) T yarr[x.size()];
  // for (std::size_t n = 0; n < x.size(); ++n) {
  //   using std::sin;
  //   yarr[n] = sin(xarr[n]);
  // }
  // const simd<T> y = loada<simd<T> >(yarr);
  // return y;
  constexpr size_t vsize = tuple_size_v<simd<T> >;
  T xarr[vsize];
  storeu(xarr, x);
  T yarr[vsize];
  for (std::size_t n = 0; n < vsize; ++n) {
    using std::sin;
    yarr[n] = sin(xarr[n]);
  }
  const simd<T> y = loadu<simd<T> >(yarr);
  return y;
}

template <typename T> struct simdl {
  typedef T value_type;
#ifndef SIMD_CPU
  nsimd::packl<T> elts;
#else
  bool elts;
#endif

  simdl(const simdl &) = default;
  simdl(simdl &&) = default;
  simdl &operator=(const simdl &) = default;
  simdl &operator=(simdl &&) = default;

  CCTK_CUDA simdl() {}
  CCTK_CUDA simdl(bool a) : elts(a) {}
#ifndef SIMD_CPU
  CCTK_CUDA simdl(const nsimd::packl<T> &elts) : elts(elts) {}
#endif

  CCTK_CUDA constexpr std::size_t size() const {
#ifndef SIMD_CPU
    return sizeof(nsimd::packl<T>) / sizeof(T);
#else
    return 1;
#endif
  }

  CCTK_CUDA simdl operator!() const { return !elts; }
  // simdl operator~() const { return ~elts; }

  friend CCTK_CUDA simdl operator&&(const simdl &x, const simdl &y) {
    return x.elts && y.elts;
  }
  friend CCTK_CUDA simdl operator&(const simdl &x, const simdl &y) {
    return x.elts && y.elts;
  }
  friend CCTK_CUDA simdl operator||(const simdl &x, const simdl &y) {
    return x.elts || y.elts;
  }
  friend CCTK_CUDA simdl operator|(const simdl &x, const simdl &y) {
    return x.elts || y.elts;
  }
  friend CCTK_CUDA simdl operator^(const simdl &x, const simdl &y) {
    return x.elts ^ y.elts;
  }

  friend CCTK_CUDA simdl operator&(const bool a, const simdl &y) {
    return a & y.elts;
  }
  friend CCTK_CUDA simdl operator|(const bool a, const simdl &y) {
    return a | y.elts;
  }
  friend CCTK_CUDA simdl operator^(const bool a, const simdl &y) {
    return a ^ y.elts;
  }

  friend CCTK_CUDA simdl operator&(const simdl &x, const bool b) {
    return x.elts & b;
  }
  friend CCTK_CUDA simdl operator|(const simdl &x, const bool b) {
    return x.elts | b;
  }
  friend CCTK_CUDA simdl operator^(const simdl &x, const bool b) {
    return x.elts ^ b;
  }

  CCTK_CUDA simdl &operator&=(const simdl &x) { return *this = *this & x; }
  CCTK_CUDA simdl &operator|=(const simdl &x) { return *this = *this | x; }
  CCTK_CUDA simdl &operator^=(const simdl &x) { return *this = *this ^ x; }

  CCTK_CUDA simdl &operator&=(const bool a) { return *this = *this & a; }
  CCTK_CUDA simdl &operator|=(const bool a) { return *this = *this | a; }
  CCTK_CUDA simdl &operator^=(const bool a) { return *this = *this ^ a; }

  friend CCTK_CUDA simdl<T> operator==(const simdl &x, const simdl &y) {
    return x.elts == y.elts;
  }
  friend CCTK_CUDA simdl<T> operator!=(const simdl &x, const simdl &y) {
    return x.elts != y.elts;
  }
  friend CCTK_CUDA simdl<T> operator<(const simdl &x, const simdl &y) {
    return x.elts < y.elts;
  }
  friend CCTK_CUDA simdl<T> operator>(const simdl &x, const simdl &y) {
    return x.elts > y.elts;
  }
  friend CCTK_CUDA simdl<T> operator<=(const simdl &x, const simdl &y) {
    return x.elts <= y.elts;
  }
  friend CCTK_CUDA simdl<T> operator>=(const simdl &x, const simdl &y) {
    return x.elts >= y.elts;
  }

  friend CCTK_CUDA simdl<T> operator==(const bool a, const simdl &y) {
    return a == y.elts;
  }
  friend CCTK_CUDA simdl<T> operator!=(const bool a, const simdl &y) {
    return a != y.elts;
  }
  friend CCTK_CUDA simdl<T> operator<(const bool a, const simdl &y) {
    return a < y.elts;
  }
  friend CCTK_CUDA simdl<T> operator>(const bool a, const simdl &y) {
    return a > y.elts;
  }
  friend CCTK_CUDA simdl<T> operator<=(const bool a, const simdl &y) {
    return a <= y.elts;
  }
  friend CCTK_CUDA simdl<T> operator>=(const bool a, const simdl &y) {
    return a >= y.elts;
  }

  friend CCTK_CUDA simdl<T> operator==(const simdl &x, const bool b) {
    return x.elts == b;
  }
  friend CCTK_CUDA simdl<T> operator!=(const simdl &x, const bool b) {
    return x.elts != b;
  }
  friend CCTK_CUDA simdl<T> operator<(const simdl &x, const bool b) {
    return x.elts < b;
  }
  friend CCTK_CUDA simdl<T> operator>(const simdl &x, const bool b) {
    return x.elts > b;
  }
  friend CCTK_CUDA simdl<T> operator<=(const simdl &x, const bool b) {
    return x.elts <= b;
  }
  friend CCTK_CUDA simdl<T> operator>=(const simdl &x, const bool b) {
    return x.elts >= b;
  }

  friend CCTK_CUDA simd<T> if_else(const simdl &cond, const simd<T> &x,
                                   const simd<T> &y) {
#ifndef SIMD_CPU
    return if_else1(cond.elts, x.elts, y.elts);
#else
    return cond.elts ? x.elts : y.elts;
#endif
  }
  friend CCTK_CUDA simd<T> if_else(const simdl &cond, const T &a,
                                   const simd<T> &y) {
#ifndef SIMD_CPU
    return if_else1(cond.elts, simd(a).elts, y.elts);
#else
    return cond.elts ? simd(a).elts : y.elts;
#endif
  }
  friend CCTK_CUDA simd<T> if_else(const simdl &cond, const simd<T> &x,
                                   const T &b) {
#ifndef SIMD_CPU
    return if_else1(cond.elts, x.elts, simd(b).elts);
#else
    return cond.elts ? x.elts : simd(b).elts;
#endif
  }
  friend CCTK_CUDA simd<T> if_else(const simdl &cond, const T &a, const T &b) {
#ifndef SIMD_CPU
    return if_else1(cond.elts, simd(a).elts, simd(b).elts);
#else
    return cond.elts ? simd(a).elts : simd(b).elts;
#endif
  }

  friend CCTK_CUDA bool all(const simdl &x) {
#ifndef SIMD_CPU
    return all(x.elts);
#else
    return x.elts;
#endif
  }
  friend CCTK_CUDA bool any(const simdl &x) {
#ifndef SIMD_CPU
    return any(x.elts);
#else
    return x.elts;
#endif
  }

  friend CCTK_CUDA void storela(T *ptr, const simdl &x) {
#ifndef SIMD_CPU
    nsimd::storela(ptr, x.elts);
#else
    *ptr = x.elts;
#endif
  }
  friend CCTK_CUDA void storelu(T *ptr, const simdl &x) {
#ifndef SIMD_CPU
    nsimd::storelu(ptr, x.elts);
#else
    *ptr = x.elts;
#endif
  }

  friend ostream &operator<<(ostream &os, const simdl &x) {
    return os << "simdl";
  }
};

} // namespace Arith
namespace std {
template <typename T>
struct tuple_size<Arith::simdl<T> >
    : std::integral_constant<std::size_t, sizeof(Arith::simdl<T>) / sizeof(T)> {
};
} // namespace std
namespace Arith {

// template <typename VT, typename T = typename VT::value_type>
// inline simdl<T> loadla(const T *ptr) {
// #ifndef SIMD_CPU
//   return nsimd::loadla<nsimd::pack<T> >(ptr);
// #else
//   return *ptr;
// #endif
// }

// template <typename VT, typename T = typename VT::value_type>
// inline simdl<T> loadlu(const T *ptr) {
// #ifndef SIMD_CPU
//   return nsimd::loadlu<nsimd::pack<T> >(ptr);
// #else
//   return *ptr;
// #endif
// }

} // namespace Arith

#endif // #ifndef SIMD_HXX
