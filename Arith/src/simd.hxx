#ifndef SIMD_HXX
#define SIMD_HXX

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
template <> struct f2i<f64> { typedef i64 type; };
template <> struct f2i<f32> { typedef i32 type; };

template <typename T> struct f2u;
template <typename T> using f2u_t = typename f2u<T>::type;
template <> struct f2u<f64> { typedef u64 type; };
template <> struct f2u<f32> { typedef u32 type; };
} // namespace detail

template <typename T> struct simd {
  typedef T value_type;
  nsimd::pack<T> elts;

  simd(const simd &) = default;
  simd(simd &&) = default;
  simd &operator=(const simd &) = default;
  simd &operator=(simd &&) = default;

  simd() {}
  simd(const T &a) : elts(a) {}
  simd(const nsimd::pack<T> &elts) : elts(elts) {}

  constexpr std::size_t size() const {
    return sizeof(nsimd::pack<T>) / sizeof(T);
  }

  simdl<T> operator!() const { return !elts; }
  simd operator~() const { return ~elts; }

  simd operator+() const { return +elts; }
  simd operator-() const { return -elts; }

  friend simd operator&(const simd &x, const simd &y) {
    return x.elts & y.elts;
  }
  friend simd operator|(const simd &x, const simd &y) {
    return x.elts | y.elts;
  }
  friend simd operator^(const simd &x, const simd &y) {
    return x.elts ^ y.elts;
  }
  friend simd operator+(const simd &x, const simd &y) {
    return x.elts + y.elts;
  }
  friend simd operator-(const simd &x, const simd &y) {
    return x.elts - y.elts;
  }
  friend simd operator*(const simd &x, const simd &y) {
    return x.elts * y.elts;
  }
  friend simd operator/(const simd &x, const simd &y) {
    return x.elts / y.elts;
  }
  friend simd operator%(const simd &x, const simd &y) {
    return x.elts % y.elts;
  }

  friend simd operator&(const T &a, const simd &y) { return a & y.elts; }
  friend simd operator|(const T &a, const simd &y) { return a | y.elts; }
  friend simd operator^(const T &a, const simd &y) { return a ^ y.elts; }
  friend simd operator+(const T &a, const simd &y) { return a + y.elts; }
  friend simd operator-(const T &a, const simd &y) { return a - y.elts; }
  friend simd operator*(const T &a, const simd &y) { return a * y.elts; }
  friend simd operator/(const T &a, const simd &y) { return a / y.elts; }
  friend simd operator%(const T &a, const simd &y) { return a % y.elts; }

  friend simd operator&(const simd &x, const T &b) { return x.elts & b; }
  friend simd operator|(const simd &x, const T &b) { return x.elts | b; }
  friend simd operator^(const simd &x, const T &b) { return x.elts ^ b; }
  friend simd operator+(const simd &x, const T &b) { return x.elts + b; }
  friend simd operator-(const simd &x, const T &b) { return x.elts - b; }
  friend simd operator*(const simd &x, const T &b) { return x.elts * b; }
  friend simd operator/(const simd &x, const T &b) { return x.elts / b; }
  friend simd operator%(const simd &x, const T &b) { return x.elts % b; }

  simd &operator&=(const simd &x) { return *this = *this & x; }
  simd &operator|=(const simd &x) { return *this = *this | x; }
  simd &operator^=(const simd &x) { return *this = *this ^ x; }
  simd &operator+=(const simd &x) { return *this = *this + x; }
  simd &operator-=(const simd &x) { return *this = *this - x; }
  simd &operator*=(const simd &x) { return *this = *this * x; }
  simd &operator/=(const simd &x) { return *this = *this / x; }
  simd &operator%=(const simd &x) { return *this = *this % x; }

  simd &operator&=(const T &a) { return *this = *this & a; }
  simd &operator|=(const T &a) { return *this = *this | a; }
  simd &operator^=(const T &a) { return *this = *this ^ a; }
  simd &operator+=(const T &a) { return *this = *this + a; }
  simd &operator-=(const T &a) { return *this = *this - a; }
  simd &operator*=(const T &a) { return *this = *this * a; }
  simd &operator/=(const T &a) { return *this = *this / a; }
  simd &operator%=(const T &a) { return *this = *this % a; }

  friend simdl<T> operator==(const simd &x, const simd &y) {
    return x.elts == y.elts;
  }
  friend simdl<T> operator!=(const simd &x, const simd &y) {
    return x.elts != y.elts;
  }
  friend simdl<T> operator<(const simd &x, const simd &y) {
    return x.elts < y.elts;
  }
  friend simdl<T> operator>(const simd &x, const simd &y) {
    return x.elts > y.elts;
  }
  friend simdl<T> operator<=(const simd &x, const simd &y) {
    return x.elts <= y.elts;
  }
  friend simdl<T> operator>=(const simd &x, const simd &y) {
    return x.elts >= y.elts;
  }

  friend simdl<T> operator==(const T &a, const simd &y) { return a == y.elts; }
  friend simdl<T> operator!=(const T &a, const simd &y) { return a != y.elts; }
  friend simdl<T> operator<(const T &a, const simd &y) { return a < y.elts; }
  friend simdl<T> operator>(const T &a, const simd &y) { return a > y.elts; }
  friend simdl<T> operator<=(const T &a, const simd &y) { return a <= y.elts; }
  friend simdl<T> operator>=(const T &a, const simd &y) { return a >= y.elts; }

  friend simdl<T> operator==(const simd &x, const T &b) { return x.elts == b; }
  friend simdl<T> operator!=(const simd &x, const T &b) { return x.elts != b; }
  friend simdl<T> operator<(const simd &x, const T &b) { return x.elts < b; }
  friend simdl<T> operator>(const simd &x, const T &b) { return x.elts > b; }
  friend simdl<T> operator<=(const simd &x, const T &b) { return x.elts <= b; }
  friend simdl<T> operator>=(const simd &x, const T &b) { return x.elts >= b; }

  friend simd abs(const simd &x) { return abs(x.elts); }
  friend simd fabs(const simd &x) { return abs(x.elts); }
  friend simdl<T> signbit(const simd &x) {
    typedef detail::f2u_t<T> U;
    constexpr U signmask = U(1) << (8 * sizeof(U) - 1);
    return to_logical(x & signmask);
  }
  friend simd sqrt(const simd &x) { return sqrt(x.elts); }
  friend simdl<T> to_logical(const simd &x) { return to_logical(x.elts); }

  friend simd copysign(const simd &x, const simd &y) {
    typedef detail::f2u_t<T> U;
    constexpr U signmask = U(1) << (8 * sizeof(U) - 1);
    return (x & ~signmask) | (y & signmask);
  }
  friend simd flipsign(const simd &x, const simd &y) {
    typedef detail::f2u_t<T> U;
    constexpr U signmask = U(1) << (8 * sizeof(U) - 1);
    return x ^ (y & signmask);
  }
  friend simd fmax(const simd &x, const simd &y) { return max(x.elts, y.elts); }
  friend simd fmin(const simd &x, const simd &y) { return min(x.elts, y.elts); }
  friend simd max(const simd &x, const simd &y) { return max(x.elts, y.elts); }
  friend simd min(const simd &x, const simd &y) { return min(x.elts, y.elts); }

  friend simd fmax(const T &a, const simd &y) { return max(a, y.elts); }
  friend simd fmin(const T &a, const simd &y) { return min(a, y.elts); }
  friend simd max(const T &a, const simd &y) { return max(a, y.elts); }
  friend simd min(const T &a, const simd &y) { return min(a, y.elts); }

  friend simd fmax(const simd &x, const T &b) { return max(x.elts, b); }
  friend simd fmin(const simd &x, const T &b) { return min(x.elts, b); }
  friend simd max(const simd &x, const T &b) { return max(x.elts, b); }
  friend simd min(const simd &x, const T &b) { return min(x.elts, b); }

  friend void storea(T *ptr, const simd &x) {
    return nsimd::storea(ptr, x.elts);
  }
  friend void storeu(T *ptr, const simd &x) {
    return nsimd::storeu(ptr, x.elts);
  }
  friend void mask_storea(const simdl<T> &mask, T *ptr, const simd &x) {
    return nsimd::mask_storea(mask.elts, ptr, x.elts);
  }
  friend void mask_storeu(const simdl<T> &mask, T *ptr, const simd &x) {
    return nsimd::mask_storeu(mask.elts, ptr, x.elts);
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
inline simd<T> iota() {
  return nsimd::iota<nsimd::pack<T> >();
}

template <typename VT, typename T = typename VT::value_type>
inline simdl<T> mask_for_loop_tail(const int i, const int n) {
  return nsimd::mask_for_loop_tail<nsimd::packl<T> >(i, n);
}

template <typename VT, typename T = typename VT::value_type>
inline simd<T> loada(const T *ptr) {
  return nsimd::loada<nsimd::pack<T> >(ptr);
}

template <typename VT, typename T = typename VT::value_type>
inline simd<T> loadu(const T *ptr) {
  return nsimd::loadu<nsimd::pack<T> >(ptr);
}

template <typename T>
inline simd<T> maskz_loada(const simdl<T> &mask, const T *ptr) {
  return nsimd::maskz_loada(mask.elts, ptr);
}

template <typename T>
inline simd<T> maskz_loadu(const simdl<T> &mask, const T *ptr) {
  return nsimd::maskz_loadu(mask.elts, ptr);
}

template <typename T>
inline simd<T> masko_loada(const simdl<T> &mask, const T *ptr,
                           const simd<T> &other) {
  return nsimd::masko_loada(mask.elts, ptr, other.elts);
}

template <typename T>
inline simd<T> masko_loadu(const simdl<T> &mask, const T *ptr,
                           const simd<T> &other) {
  return nsimd::masko_loadu(mask.elts, ptr, other.elts);
}

template <typename T, typename U,
          enable_if_t<is_convertible_v<T, U> > * = nullptr>
inline simd<T> masko_loada(const simdl<T> &mask, const T *ptr, const U &other) {
  return masko_loada(mask, ptr, simd<T>(other));
}

template <typename T, typename U,
          enable_if_t<is_convertible_v<T, U> > * = nullptr>
inline simd<T> masko_loadu(const simdl<T> &mask, const T *ptr, const U &other) {
  return masko_loadu(mask, ptr, simd<T>(other));
}

template <typename T> inline simd<T> cbrt(const simd<T> &x) {
  // alignas(alignof(simd<T>)) T xarr[x.size()];
  // storea(xarr, x);
  // alignas(alignof(simd<T>)) T yarr[x.size()];
  // for (std::size_t n = 0; n < x.size(); ++n) {
  //   using std::cbrt;
  //   yarr[n] = cbrt(xarr[n]);
  // }
  // const simd<T> y = loada<simd<T> >(yarr);
  // return y;
  T xarr[x.size()];
  storeu(xarr, x);
  T yarr[x.size()];
  for (std::size_t n = 0; n < x.size(); ++n) {
    using std::cbrt;
    yarr[n] = cbrt(xarr[n]);
  }
  const simd<T> y = loadu<simd<T> >(yarr);
  return y;
}

template <typename T> inline simd<T> sin(const simd<T> &x) {
  // alignas(alignof(simd<T>)) T xarr[x.size()];
  // storea(xarr, x);
  // alignas(alignof(simd<T>)) T yarr[x.size()];
  // for (std::size_t n = 0; n < x.size(); ++n) {
  //   using std::sin;
  //   yarr[n] = sin(xarr[n]);
  // }
  // const simd<T> y = loada<simd<T> >(yarr);
  // return y;
  T xarr[x.size()];
  storeu(xarr, x);
  T yarr[x.size()];
  for (std::size_t n = 0; n < x.size(); ++n) {
    using std::sin;
    yarr[n] = sin(xarr[n]);
  }
  const simd<T> y = loadu<simd<T> >(yarr);
  return y;
}

template <typename T> struct simdl {
  typedef T value_type;
  nsimd::packl<T> elts;

  simdl(const simdl &) = default;
  simdl(simdl &&) = default;
  simdl &operator=(const simdl &) = default;
  simdl &operator=(simdl &&) = default;

  simdl() {}
  simdl(T a) : elts(a) {}
  simdl(const nsimd::packl<T> &elts) : elts(elts) {}

  constexpr std::size_t size() const {
    return sizeof(nsimd::packl<T>) / sizeof(T);
  }

  simdl operator!() const { return !elts; }
  simdl operator~() const { return ~elts; }

  friend simdl operator&&(const simdl &x, const simdl &y) {
    return x.elts && y.elts;
  }
  friend simdl operator&(const simdl &x, const simdl &y) {
    return x.elts && y.elts;
  }
  friend simdl operator||(const simdl &x, const simdl &y) {
    return x.elts || y.elts;
  }
  friend simdl operator|(const simdl &x, const simdl &y) {
    return x.elts || y.elts;
  }
  friend simdl operator^(const simdl &x, const simdl &y) {
    return x.elts ^ y.elts;
  }

  friend simdl operator&(const bool a, const simdl &y) { return a & y.elts; }
  friend simdl operator|(const bool a, const simdl &y) { return a | y.elts; }
  friend simdl operator^(const bool a, const simdl &y) { return a ^ y.elts; }

  friend simdl operator&(const simdl &x, const bool b) { return x.elts & b; }
  friend simdl operator|(const simdl &x, const bool b) { return x.elts | b; }
  friend simdl operator^(const simdl &x, const bool b) { return x.elts ^ b; }

  simdl &operator&=(const simdl &x) { return *this = *this & x; }
  simdl &operator|=(const simdl &x) { return *this = *this | x; }
  simdl &operator^=(const simdl &x) { return *this = *this ^ x; }

  simdl &operator&=(const bool a) { return *this = *this & a; }
  simdl &operator|=(const bool a) { return *this = *this | a; }
  simdl &operator^=(const bool a) { return *this = *this ^ a; }

  friend simdl<T> operator==(const simdl &x, const simdl &y) {
    return x.elts == y.elts;
  }
  friend simdl<T> operator!=(const simdl &x, const simdl &y) {
    return x.elts != y.elts;
  }
  friend simdl<T> operator<(const simdl &x, const simdl &y) {
    return x.elts < y.elts;
  }
  friend simdl<T> operator>(const simdl &x, const simdl &y) {
    return x.elts > y.elts;
  }
  friend simdl<T> operator<=(const simdl &x, const simdl &y) {
    return x.elts <= y.elts;
  }
  friend simdl<T> operator>=(const simdl &x, const simdl &y) {
    return x.elts >= y.elts;
  }

  friend simdl<T> operator==(const bool a, const simdl &y) {
    return a == y.elts;
  }
  friend simdl<T> operator!=(const bool a, const simdl &y) {
    return a != y.elts;
  }
  friend simdl<T> operator<(const bool a, const simdl &y) { return a < y.elts; }
  friend simdl<T> operator>(const bool a, const simdl &y) { return a > y.elts; }
  friend simdl<T> operator<=(const bool a, const simdl &y) {
    return a <= y.elts;
  }
  friend simdl<T> operator>=(const bool a, const simdl &y) {
    return a >= y.elts;
  }

  friend simdl<T> operator==(const simdl &x, const bool b) {
    return x.elts == b;
  }
  friend simdl<T> operator!=(const simdl &x, const bool b) {
    return x.elts != b;
  }
  friend simdl<T> operator<(const simdl &x, const bool b) { return x.elts < b; }
  friend simdl<T> operator>(const simdl &x, const bool b) { return x.elts > b; }
  friend simdl<T> operator<=(const simdl &x, const bool b) {
    return x.elts <= b;
  }
  friend simdl<T> operator>=(const simdl &x, const bool b) {
    return x.elts >= b;
  }

  friend simd<T> if_else(const simdl &cond, const simd<T> &x,
                         const simd<T> &y) {
    return if_else1(cond.elts, x.elts, y.elts);
  }
  friend simd<T> if_else(const simdl &cond, const T &a, const simd<T> &y) {
    return if_else1(cond.elts, simd(a).elts, y.elts);
  }
  friend simd<T> if_else(const simdl &cond, const simd<T> &x, const T &b) {
    return if_else1(cond.elts, x.elts, simd(b).elts);
  }
  friend simd<T> if_else(const simdl &cond, const T &a, const T &b) {
    return if_else1(cond.elts, simd(a).elts, simd(b).elts);
  }

  friend bool all(const simdl &x) { return all(x.elts); }
  friend bool any(const simdl &x) { return any(x.elts); }

  friend void storela(T *ptr, const simdl &x) {
    return nsimd::storela(ptr, x.elts);
  }
  friend void storelu(T *ptr, const simdl &x) {
    return nsimd::storelu(ptr, x.elts);
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

template <typename VT, typename T = typename VT::value_type>
inline simdl<T> loadla(const T *ptr) {
  return nsimd::loadla<nsimd::pack<T> >(ptr);
}

template <typename VT, typename T = typename VT::value_type>
inline simdl<T> loadlu(const T *ptr) {
  return nsimd::loadlu<nsimd::pack<T> >(ptr);
}

} // namespace Arith

#endif // #ifndef SIMD_HXX
