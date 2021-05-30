#ifndef DEFS_HXX
#define DEFS_HXX

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>

#include <cmath>
#include <limits>
#include <type_traits>

#ifdef CCTK_DEBUG
#define ARITH_INLINE
#else
#define ARITH_INLINE CCTK_ATTRIBUTE_ALWAYS_INLINE
#endif

#ifdef __CUDACC__
#define ARITH_DEVICE __device__
#define ARITH_HOST __host__
#else
#define ARITH_DEVICE
#define ARITH_HOST
#endif

namespace Arith {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// Return the value zero for a given type
template <typename T> struct zero;
// template <typename T> inline constexpr T zero_v = zero<T>::value;

template <> struct zero<bool> : integral_constant<bool, 0> {};
template <> struct zero<short> : integral_constant<short, 0> {};
template <> struct zero<int> : integral_constant<int, 0> {};
template <> struct zero<long> : integral_constant<long, 0> {};
template <> struct zero<long long> : integral_constant<long long, 0> {};

template <> struct zero<float> {
  typedef float value_type;
  static constexpr value_type value = 0;
  constexpr ARITH_INLINE operator value_type() const { return value; }
  constexpr ARITH_INLINE value_type operator()() const { return value; }
};

template <> struct zero<double> {
  typedef double value_type;
  static constexpr value_type value = 0;
  constexpr ARITH_INLINE operator value_type() const { return value; }
  constexpr ARITH_INLINE value_type operator()() const { return value; }
};

// Return the value one for a given type
template <typename T> struct one;
// template <typename T> inline constexpr T one_v = one<T>::value;

template <> struct one<bool> : integral_constant<bool, 1> {};
template <> struct one<short> : integral_constant<short, 1> {};
template <> struct one<int> : integral_constant<int, 1> {};
template <> struct one<long> : integral_constant<long, 1> {};

template <> struct one<float> {
  typedef float value_type;
  static constexpr value_type value = 1;
  constexpr ARITH_INLINE operator value_type() const { return value; }
  constexpr ARITH_INLINE value_type operator()() const { return value; }
};

template <> struct one<double> {
  typedef double value_type;
  static constexpr value_type value = 1;
  constexpr ARITH_INLINE operator value_type() const { return value; }
  constexpr ARITH_INLINE value_type operator()() const { return value; }
};

// Return the value nan for a given type
template <typename T> struct nan;
// template <typename T> inline constexpr T nan_v = nan<T>::value;

template <>
struct nan<bool> : integral_constant<bool, numeric_limits<bool>::max()> {};
template <>
struct nan<short> : integral_constant<short, numeric_limits<short>::min()> {};
template <>
struct nan<int> : integral_constant<int, numeric_limits<int>::min()> {};
template <>
struct nan<long> : integral_constant<long, numeric_limits<long>::min()> {};

template <> struct nan<float> {
  typedef float value_type;
  static constexpr value_type value = numeric_limits<value_type>::quiet_NaN();
  constexpr ARITH_INLINE operator value_type() const { return value; }
  constexpr ARITH_INLINE value_type operator()() const { return value; }
};

template <> struct nan<double> {
  typedef double value_type;
  static constexpr value_type value = numeric_limits<value_type>::quiet_NaN();
  constexpr ARITH_INLINE operator value_type() const { return value; }
  constexpr ARITH_INLINE value_type operator()() const { return value; }
};

////////////////////////////////////////////////////////////////////////////////

template <int imin, int imax, int istep = 1, typename F,
          enable_if_t<(istep > 0 ? imin >= imax : imin <= imax)> * = nullptr>
ARITH_INLINE void unroll_for(const F &f) {
  // done: do nothing
}

template <int imin, int imax, int istep = 1, typename F,
          enable_if_t<!(istep > 0 ? imin >= imax : imin <= imax)> * = nullptr>
ARITH_INLINE void unroll_for(const F &f) {
  f(imin);
  unroll_for<imin + istep, imax, istep>(f);
}

////////////////////////////////////////////////////////////////////////////////

// Return true if all elements of a container are true
constexpr ARITH_INLINE bool all(bool x) { return x; }

// Return true if any elements of a container are true
constexpr ARITH_INLINE bool any(bool x) { return x; }

// Return true if there is a nan anywhere in x. This always returns a bool, even
// when x is a vector or SIMD type that would otherwise perform operations
// pointwise. Use this for debugging.
constexpr ARITH_INLINE bool anyisnan(const float &x) {
  using std::isnan;
  return isnan(x);
}
constexpr ARITH_INLINE bool anyisnan(const double &x) {
  using std::isnan;
  return isnan(x);
}

// if-then-else function that can be used with SIMD types
template <typename T, typename U>
constexpr ARITH_INLINE auto if_else(bool c, const T &x, const U &y) {
  return c ? x : y;
}

////////////////////////////////////////////////////////////////////////////////

// bitsign(i) = (-1)^i, the converse to signbit
constexpr int bitsign(bool c) { return if_else(c, -1, 1); }
constexpr int bitsign(int i) { return bitsign(i % 2 != 0); }

// Return x if y>0, -x if y<0
template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T flipsign(const T &x, const T &y) {
  using std::copysign;
  return copysign(T(1), y) * x;
}

// A max function that returns nan when any argument is nan
template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T max1(const T &x, const T &y) {
  using std::max;
  return if_else(x != x, x, if_else(y != y, y, max(x, y)));
}

// The maximum of the absolute values. This is reduces over containers.
template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T maxabs(const T &x) {
  using std::abs;
  return abs(x);
}

// A min function that returns nan when any argument is nan
template <typename T>
inline ARITH_INLINE ARITH_DEVICE ARITH_HOST T min1(const T &x, const T &y) {
  using std::min;
  return if_else(x != x, x, if_else(y != y, y, min(x, y)));
}

namespace detail {
template <typename T> constexpr T pown(const T &x, int n) {
  T r{1};
  T y{x};
  // invariant: initial(x^n) == r * y^n
  while (n) {
    if (n & 1)
      r *= y;
    y *= y;
    n >>= 1;
  }
  return r;
}
} // namespace detail

// Raise a value to an integer power, calculated efficiently
template <typename T> constexpr T pown(const T &x, const int n) {
  return n >= 0 ? detail::pown(x, n) : 1 / detail::pown(x, -n);
}

template <typename T> constexpr T pow2(const T &x) { return x * x; }

} // namespace Arith

#endif // #ifndef DEFS_HXX
