#ifndef SIMT_HXX
#define SIMT_HXX

// We need this before including nsimd
#include <cctk.h>

#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#include <nsimd/nsimd-all.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>

namespace SIMT {

////////////////////////////////////////////////////////////////////////////////

// A few additional standard functions that are useful when using vectors

using std::abs;
using std::max;
using std::min;

template <typename T, typename F>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, void>
forall(std::ptrdiff_t imin, std::ptrdiff_t imax, const F &f) {
  for (std::ptrdiff_t i = imin; i < imax; ++i)
    f(i, true);
}

template <typename T>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> iota(const T &) {
  return T(0);
}

template <typename T>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, int> len(const T &) {
  return 1;
}
template <typename T>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, int> len() {
  return 1;
}

template <typename T>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> maxval(const T &x) {
  return x;
}

template <typename T>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> minval(const T &x) {
  return x;
}

template <typename T>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> rec(const T &x) {
  return T(1) / x;
}

template <typename T, typename Op, typename F>
constexpr T reduce(const Op &op, T res, std::ptrdiff_t imin,
                   std::ptrdiff_t imax, std::ptrdiff_t istep, const F &f) {
  for (std::ptrdiff_t i = imin; i < imax; i += istep)
    res = op(res, f(i));
  return res;
}
template <typename T, typename Op, typename F>
constexpr T reduce(const Op &op, T res, std::ptrdiff_t imin,
                   std::ptrdiff_t imax, const F &f) {
  return reduce(op, res, imin, imax, 1, f);
}
template <typename T, typename Op, typename F>
constexpr T reduce(const Op &op, T res, std::ptrdiff_t imax, const F &f) {
  return reduce(op, res, 0, imax, f);
}

template <typename T>
constexpr std::enable_if_t<std::is_arithmetic_v<T>, T> sum(const T &x) {
  return x;
}

////////////////////////////////////////////////////////////////////////////////

// A reference or a pointer to memory

template <typename T, typename U> struct memref_assign;

template <typename T> struct memref {
  T &ref;
  memref() = delete;
  constexpr memref(T &ref) : ref(ref) {}
  template <typename U> constexpr memref &operator=(const U &x) {
    memref_assign<T, U>()(*this, x);
    return *this;
  }

  template <typename U> constexpr memref &operator+=(const U &x) {
    return *this = *this + x;
  }
  template <typename U> constexpr memref &operator-=(const U &x) {
    return *this = *this - x;
  }
  template <typename U> constexpr memref &operator*=(const U &x) {
    return *this = *this * x;
  }
  template <typename U> constexpr memref &operator/=(const U &x) {
    return *this = *this / x;
  }
};

template <typename T> constexpr memref<T> make_memref(T &ref) { return ref; }

template <typename T> struct memptr {
  T *ptr;
  memptr() = default;
  constexpr memptr(T *ptr) : ptr(ptr) {}
  constexpr memptr(std::intptr_t ptr) : ptr(ptr) {}
  constexpr memptr(std::uintptr_t ptr) : ptr(ptr) {}
  constexpr operator T *() const { return ptr; }
  constexpr operator std::intptr_t() const { return ptr; }
  constexpr operator std::uintptr_t() const { return ptr; }

  constexpr memref<T> operator*() const { return *ptr; }
  constexpr memref<T> operator[](std::ptrdiff_t i) const { return ptr[i]; }
  constexpr memref<T> operator->() const { return ptr; }

  constexpr operator bool() const { return ptr; }
  constexpr memptr &operator++() { return ++ptr, *this; }
  constexpr memptr operator++(int) { return ptr++; }
  constexpr memptr &operator--() { return --ptr, *this; }
  constexpr memptr operator--(int) { return ptr--; }
  friend constexpr memptr operator+(const memptr &mptr,
                                    const std::ptrdiff_t &i) {
    return mptr.ptr + i;
  }
  friend constexpr memptr operator+(const std::ptrdiff_t &i,
                                    const memptr &mptr) {
    return i + mptr.ptr;
  }
  friend constexpr memptr operator-(const memptr &mptr,
                                    const std::ptrdiff_t &i) {
    return mptr.ptr + i;
  }
  friend constexpr std::ptrdiff_t operator-(const memptr &mptr1,
                                            const memptr &mptr2) {
    return mptr1.ptr - mptr2.ptr;
  }
};

template <typename T> constexpr memptr<T> make_memptr(T *ptr) { return ptr; }

////////////////////////////////////////////////////////////////////////////////

// The simt class

template <typename T> struct simt {
  typedef nsimd::pack<T> VT;
  // static constexpr size_t N = nsimd::len<VT>();
  static constexpr size_t N = sizeof(VT) / sizeof(T);

  VT elems;

  typedef T value_type;
  constexpr size_t size() const { return N; }

  simt() = default;

  constexpr simt(const VT &a) : elems(a) {}
  constexpr simt(const std::array<T, N> &arr)
      : elems(nsimd::loadu<VT>(arr.data())) {}
  constexpr std::array<T, N> to_array() const {
    std::array<T, N> arr;
    nsimd::storeu(arr.data(), elems);
    return arr;
  }

  constexpr simt(const T &a) : elems(nsimd::set1<VT>(a)) {}

  constexpr simt(const memref<T> &mref) : elems(nsimd::loadu<VT>(&mref.ref)) {}

  friend constexpr simt operator+(const simt &x) { return simt(x.elems); }
  friend constexpr simt operator-(const simt &x) { return simt(-x.elems); }

  friend constexpr simt operator+(const simt &x, const simt &y) {
    return simt(x.elems + y.elems);
  }
  friend constexpr simt operator-(const simt &x, const simt &y) {
    return simt(x.elems - y.elems);
  }
  friend constexpr simt operator*(const simt &x, const simt &y) {
    return simt(x.elems * y.elems);
  }
  friend constexpr simt operator/(const simt &x, const simt &y) {
    return simt(x.elems / y.elems);
  }

  constexpr simt &operator+=(const simt &x) { return *this = *this + x; }
  constexpr simt &operator-=(const simt &x) { return *this = *this - x; }
  constexpr simt &operator*=(const simt &x) { return *this = *this * x; }
  constexpr simt &operator/=(const simt &x) { return *this = *this / x; }

  friend std::ostream &operator<<(std::ostream &os, const simt &x) {
    return os << x.elems;
  }
};

////////////////////////////////////////////////////////////////////////////////

// Stand-alone simt functions

// Declarations

template <typename T> constexpr simt<T> abs(const simt<T> &x);
template <typename T> constexpr simt<T> fabs(const simt<T> &x);

template <typename T, typename F>
constexpr
    std::enable_if_t<std::is_same_v<T, simt<typename T::value_type> >, void>
    forall(std::ptrdiff_t imin, std::ptrdiff_t imax, const F &f);

// TODO: Introduce template specialization for simt<bool>
template <typename T>
constexpr T ifelse(const simt<bool> &cond, const simt<T> &x, const simt<T> &y);

template <typename T> constexpr simt<T> iota(const simt<T> &);

template <typename T> constexpr int len(const simt<T> &);
template <typename T>
constexpr
    std::enable_if_t<std::is_same_v<T, simt<typename T::value_type> >, int>
    len();

template <typename T> constexpr simt<T> max(const simt<T> &x, const simt<T> &y);
template <typename T>
constexpr simt<T> fmax(const simt<T> &x, const simt<T> &y);

template <typename T> constexpr T maxval(const simt<T> &x);

template <typename T> constexpr simt<T> min(const simt<T> &x, const simt<T> &y);
template <typename T>
constexpr simt<T> fmin(const simt<T> &x, const simt<T> &y);

template <typename T> constexpr T minval(const simt<T> &x);

template <typename T, typename Gen, typename Dist>
simt<T> random(const simt<T> &, Gen &gen, Dist &dist);

template <typename T> constexpr simt<T> rec(const simt<T> &x);

template <typename T> constexpr T sum(const simt<T> &x);

// Definitions

template <typename T> constexpr simt<T> abs(const simt<T> &x) {
  return simt<T>(abs(x.elems));
}
template <typename T> constexpr simt<T> fabs(const simt<T> &x) {
  return abs(x);
}

template <typename T, typename F>
constexpr
    std::enable_if_t<std::is_same_v<T, simt<typename T::value_type> >, void>
    forall(std::ptrdiff_t imin, std::ptrdiff_t imax, const F &f) {
  // typedef nsimd::pack<CCTK_INT8> I;
  // typedef nsimd::pack<long int> I;
  // typedef nsimd::pack<long long int> I;
  // typedef nsimd::pack<nsimd::longlong> I;
  typedef nsimd::pack<i64> I;
  // typedef nsimd::packl<typename T::value_type> B;
  typedef nsimd::packl<typename I::value_type> B;
  constexpr int vlen = len<T>();
  for (std::ptrdiff_t i = imin; i < imax; i += vlen) {
    B mask =
        nsimd::iota<I>() < nsimd::set1<I>((typename I::value_type)(imax - i));
    nsimd::packl<typename T::value_type> mask2 =
        nsimd::reinterpretl<nsimd::packl<typename T::value_type> >(mask);
    f(i, mask);
  }
}

template <typename T>
constexpr T ifelse(const simt<bool> &cond, const simt<T> &x, const simt<T> &y) {
  return nsimd::if_else(cond, x, y);
}

template <typename T> constexpr simt<T> iota(const simt<T> &) {
  return simt<T>(nsimd::iota<nsimd::pack<T> >());
}

template <typename T> constexpr int len(const simt<T> &) { return simt<T>::N; }
template <typename T>
constexpr
    std::enable_if_t<std::is_same_v<T, simt<typename T::value_type> >, int>
    len() {
  return T::N;
}

template <typename T>
constexpr simt<T> max(const simt<T> &x, const simt<T> &y) {
  return simt<T>(max(x.elems, y.elems));
}
template <typename T>
constexpr simt<T> fmax(const simt<T> &x, const simt<T> &y) {
  return max(x, y);
}

template <typename T> constexpr T maxval(const simt<T> &x) {
  // const auto arr = x.to_array();
  // return *std::max_element(arr.begin(), arr.end());
  simt<T> dummy;
  constexpr std::size_t N = dummy.size();
  simt<T> r = x;
  for (std::size_t n = N; n > 1; n /= 2) {
    simt<T> lo = nsimd::unziplo(r.elems, r.elems);
    simt<T> hi = nsimd::unziphi(r.elems, r.elems);
    r = max(lo, hi);
  }
  return r.to_array()[0];
}

template <typename T>
constexpr simt<T> min(const simt<T> &x, const simt<T> &y) {
  return simt<T>(min(x.elems, y.elems));
}
template <typename T>
constexpr simt<T> fmin(const simt<T> &x, const simt<T> &y) {
  return min(x, y);
}

template <typename T> constexpr T minval(const simt<T> &x) {
  // const auto arr = x.to_array();
  // return *std::min_element(arr.begin(), arr.end());
  simt<T> dummy;
  constexpr std::size_t N = dummy.size();
  simt<T> r = x;
  for (std::size_t n = N; n > 1; n /= 2) {
    simt<T> lo = nsimd::unziplo(r.elems, r.elems);
    simt<T> hi = nsimd::unziphi(r.elems, r.elems);
    r = min(lo, hi);
  }
  return r.to_array()[0];
}

template <typename T, typename Gen, typename Dist>
simt<T> random(const simt<T> &, Gen &gen, Dist &dist) {
  SIMT::simt<T> dummy;
  constexpr std::size_t N = dummy.size();
  std::array<T, N> arr;
  for (T &elt : arr)
    elt = dist(gen);
  return simt(arr);
}

template <typename T> constexpr simt<T> rec(const simt<T> &x) {
  return simt<T>(rec(x.elems));
}

template <typename T> constexpr T sum(const simt<T> &x) {
  return addv(x.elems);
}

////////////////////////////////////////////////////////////////////////////////

// Template specializations for SIMT::simt

template <typename T> struct memref_assign<T, simt<T> > {
  constexpr memref<T> operator()(memref<T> mref, const simt<T> &x) const {
    nsimd::storeu(&mref.ref, x.elems);
    return mref;
  }
};

} // namespace SIMT

namespace std {

template <typename T> struct equal_to<SIMT::simt<T> > {
  constexpr bool operator()(const SIMT::simt<T> &lhs,
                            const SIMT::simt<T> &rhs) const {
    constexpr std::size_t N = std::tuple_size_v<SIMT::simt<T> >;
    return equal_to<std::array<T, N> >()(lhs.to_array(), rhs.to_array());
  }
};

template <std::size_t I, typename T> constexpr T get(const SIMT::simt<T> &x) {
  constexpr std::size_t N = SIMT::simt<T>::N;
  return tuple_element<I, std::array<T, N> >()(x.to_array());
}

template <typename T> class numeric_limits<SIMT::simt<T> > {
public:
  static constexpr SIMT::simt<T> min() { return numeric_limits<T>::min(); }
  static constexpr SIMT::simt<T> lowest() {
    return numeric_limits<T>::lowest();
  }
  static constexpr SIMT::simt<T> max() { return numeric_limits<T>::max(); }
  static constexpr SIMT::simt<T> epsilon() {
    return numeric_limits<T>::epsilon();
  }
  static constexpr SIMT::simt<T> round_error() {
    return numeric_limits<T>::round_error();
  }
  static constexpr SIMT::simt<T> quiet_NaN() {
    return numeric_limits<T>::quiet_NaN();
  }
  static constexpr SIMT::simt<T> signaling_NaN() {
    return numeric_limits<T>::signaling_NaN();
  }
  static constexpr SIMT::simt<T> denorm_min() {
    return numeric_limits<T>::denorm_min();
  }
};

template <std::size_t I, typename T> struct tuple_element<I, SIMT::simt<T> > {
  typedef T type;
};

template <typename T>
struct tuple_size<SIMT::simt<T> >
    : std::integral_constant<std::size_t, SIMT::simt<T>::N> {};

} // namespace std

#endif // #ifndef SIMT_HXX
