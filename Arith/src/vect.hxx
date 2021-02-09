#ifndef VECT_HXX
#define VECT_HXX

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>

#include <array>
#include <complex>
#include <initializer_list>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace Arith {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

template <typename T, size_t N> struct ntuple {
  typedef decltype(tuple_cat(declval<tuple<T> >(),
                             declval<typename ntuple<T, N - 1>::type>())) type;
};
template <typename T> struct ntuple<T, 0> { typedef tuple<> type; };
template <typename T, size_t N> using ntuple_t = typename ntuple<T, N>::type;

////////////////////////////////////////////////////////////////////////////////

// TODO: This is broken (at least on MacOS, with gcc 9.2, with -O0).

// namespace {
// template <typename T, size_t N, typename... Ts>
// CCTK_ATTRIBUTE_ALWAYS_INLINE
// constexpr enable_if_t<(sizeof...(Ts) == N), array<T, N> >
// array_from_initializer_list(const T *const beg, const T *const end,
//                             const Ts &... xs) {
//   return array<T, N>{xs...};
// }
//
// template <typename T, size_t N, typename... Ts>
// CCTK_ATTRIBUTE_ALWAYS_INLINE
// constexpr enable_if_t<(sizeof...(Ts) < N), array<T, N> >
// array_from_initializer_list(const T *const beg, const T *const end,
//                             const Ts &... xs) {
//   return array_from_initializer_list<T, N>(beg + 1, end, xs...,
//                                            beg < end ? *beg : T{});
// }
// } // namespace
//
// template <typename T, size_t N>
//  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE array<T, N>
// array_from_initializer_list(initializer_list<T> l) {
//   return array_from_initializer_list<T, N>(l.begin(), l.end());
// }

// namespace {
// template <typename T, size_t N, typename... Ts>
// CCTK_ATTRIBUTE_ALWAYS_INLINE
// constexpr enable_if_t<(sizeof...(Ts) == N), array<T, N> >
// array_from_initializer_list(const T *const beg, const T *const end,
//                             const Ts &... xs) {
//   return array<T, N>{xs...};
// }

// template <typename T, size_t N, typename... Ts>
// CCTK_ATTRIBUTE_ALWAYS_INLINE
// constexpr enable_if_t<(sizeof...(Ts) < N), array<T, N> >
// array_from_initializer_list(const T *const beg, const T *const end,
//                             const Ts &... xs) {
//   return array_from_initializer_list<T, N>(beg, end - 1,
//                                            beg < end ? *(end - 1) : T{},
//                                            xs...);
// }
// } // namespace
//
// template <typename T, size_t N>
//  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE array<T, N>
// array_from_initializer_list(initializer_list<T> l) {
//   return array_from_initializer_list<T, N>(l.begin(), l.end());
// }

template <typename T, size_t N>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE enable_if_t<(N == 3), array<T, N> >
array_from_initializer_list(initializer_list<T> l) {
#ifndef __CUDACC__
  assert(l.size() >= N);
#endif
  const T *restrict const p = l.begin();
  return {p[0], p[1], p[2]};
}

template <typename T, size_t N>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE enable_if_t<(N == 4), array<T, N> >
array_from_initializer_list(initializer_list<T> l) {
#ifndef __CUDACC__
  assert(l.size() >= N);
#endif
  const T *restrict const p = l.begin();
  return {p[0], p[1], p[2], p[3]};
}

template <typename T, size_t N>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE enable_if_t<(N == 6), array<T, N> >
array_from_initializer_list(initializer_list<T> l) {
#ifndef __CUDACC__
  assert(l.size() >= N);
#endif
  const T *restrict const p = l.begin();
  return {p[0], p[1], p[2], p[3], p[4], p[5]};
}

template <typename T, size_t N>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE enable_if_t<(N == 10), array<T, N> >
array_from_initializer_list(initializer_list<T> l) {
#ifndef __CUDACC__
  assert(l.size() >= N);
#endif
  const T *restrict const p = l.begin();
  return {p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], p[9]};
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE array<T, 3>
array_from_tuple(tuple<T, T, T> t) {
  return {move(get<0>(t)), move(get<1>(t)), move(get<2>(t))};
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE array<T, 4>
array_from_tuple(tuple<T, T, T, T> t) {
  return {move(get<0>(t)), move(get<1>(t)), move(get<2>(t)), move(get<3>(t))};
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE array<T, 6>
array_from_tuple(tuple<T, T, T, T, T, T> t) {
  return {move(get<0>(t)), move(get<1>(t)), move(get<2>(t)),
          move(get<3>(t)), move(get<4>(t)), move(get<5>(t))};
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE array<T, 10>
array_from_tuple(tuple<T, T, T, T, T, T, T, T, T, T> t) {
  return {move(get<0>(t)), move(get<1>(t)), move(get<2>(t)), move(get<3>(t)),
          move(get<4>(t)), move(get<5>(t)), move(get<6>(t)), move(get<7>(t)),
          move(get<8>(t)), move(get<9>(t))};
}

template <typename T, size_t N>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE enable_if_t<(N == 3), array<T, N> >
array_from_vector(vector<T> &&v) {
#ifndef __CUDACC__
  assert(v.size() >= N);
#endif
  typename vector<T>::iterator const p = v.begin();
  return {move(p[0]), move(p[1]), move(p[2])};
}

template <typename T, size_t N>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE enable_if_t<(N == 4), array<T, N> >
array_from_vector(vector<T> &&v) {
#ifndef __CUDACC__
  assert(v.size() >= N);
#endif
  typename vector<T>::iterator const p = v.begin();
  return {move(p[0]), move(p[1]), move(p[2]), move(p[3])};
}

template <typename T, size_t N>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE enable_if_t<(N == 6), array<T, N> >
array_from_vector(vector<T> &&v) {
#ifndef __CUDACC__
  assert(v.size() >= N);
#endif
  typename vector<T>::iterator const p = v.begin();
  return {move(p[0]), move(p[1]), move(p[2]),
          move(p[3]), move(p[4]), move(p[5])};
}

template <typename T, size_t N>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE enable_if_t<(N == 10), array<T, N> >
array_from_vector(vector<T> &&v) {
#ifndef __CUDACC__
  assert(v.size() >= N);
#endif
  typename vector<T>::iterator const p = v.begin();
  return {move(p[0]), move(p[1]), move(p[2]), move(p[3]), move(p[4]),
          move(p[5]), move(p[6]), move(p[7]), move(p[8]), move(p[9])};
}

static_assert(static_cast<const array<int, 3> &>(
                  array_from_initializer_list<int, 3>({0, 1, 2}))[0] == 0,
              "");
static_assert(static_cast<const array<int, 3> &>(
                  array_from_initializer_list<int, 3>({0, 1, 2}))[1] == 1,
              "");
static_assert(static_cast<const array<int, 3> &>(
                  array_from_initializer_list<int, 3>({0, 1, 2}))[2] == 2,
              "");

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct zero;

template <> struct zero<bool> {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE bool operator()() const { return 0; }
};

template <> struct zero<short> {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE short operator()() const { return 0; }
};

template <> struct zero<int> {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE int operator()() const { return 0; }
};

template <> struct zero<long> {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE long operator()() const { return 0; }
};

template <> struct zero<float> {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE float operator()() const { return 0; }
};

template <> struct zero<double> {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE double operator()() const { return 0; }
};

template <typename T> struct zero<complex<T> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE complex<T> operator()() const {
    return complex<T>(zero<T>()(), zero<T>()());
  }
};

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D> struct vect {
  array<T, D> elts;

  // initializes all elements to zero
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect() : elts() {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(const vect &) = default;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(vect &&) = default;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect &
  operator=(const vect &) = default;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect &operator=(vect &&) = default;

  template <typename U>
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(vect<U, D> x) : elts() {
    for (int d = 0; d < D; ++d)
      elts[d] = move(x.elts[d]);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(array<T, D> arr)
      : elts(move(arr)) {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(ntuple_t<T, D> tup)
      : elts(array_from_tuple(move(tup))) {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(initializer_list<T> lst)
      : elts(array_from_initializer_list<T, D>(lst)) {
#ifdef CCTK_DEBUG
    assert(lst.size() == D);
#endif
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(const T (&arr)[D]) : elts() {
    for (int d = 0; d < D; ++d)
      elts[d] = arr[d];
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(const vector<T> &vec)
      : elts(array_from_vector<T, D>(vec)) {
#ifdef CCTK_DEBUG
    assert(vec.size() == D);
#endif
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect(vector<T> &&vec)
      : elts(array_from_vector<T, D>(move(vec))) {}

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect pure(T a) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = a;
    return r;
  }

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect unit(int dir) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = d == dir;
    return r;
  }

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect iota() {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = d;
    return r;
  }

  template <typename F>
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect make(const F &f) {
    return vect<int, D>::iota().map(f);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE const T &operator[](int d) const {
#ifdef CCTK_DEBUG
    assert(d >= 0 && d < D);
#endif
    return elts[d];
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE T &operator[](int d) {
#ifdef CCTK_DEBUG
    assert(d >= 0 && d < D);
#endif
    return elts[d];
  }

  template <typename F,
            typename R = remove_cv_t<remove_reference_t<result_of_t<F(T)> > > >
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<R, D> map(const F &f) const {
    vect<R, D> r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = f(elts[d]);
    return r;
  }
  template <
      typename F, typename U,
      typename R = remove_cv_t<remove_reference_t<result_of_t<F(T, U)> > > >
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<R, D>
  map2(const F &f, const vect<U, D> &y) const {
    vect<R, D> r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = f(elts[d], y.elts[d]);
    return r;
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator+(const vect &x) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = +x.elts[d];
    return r;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator-(const vect &x) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = -x.elts[d];
    return r;
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator+(const vect &x,
                                                               const vect &y) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] + y.elts[d];
    return r;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator-(const vect &x,
                                                               const vect &y) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] - y.elts[d];
    return r;
  }

  // friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect
  // operator+(const vect &x, const array<T, D> &y) {
  //   return x + vect(y);
  // }
  // friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect
  // operator-(const vect &x, const array<T, D> &y) {
  //   return x - vect(y);
  // }

  // friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator+(const vect &x,
  //                                                              T a) {
  //   vect r;
  //   for (int d = 0; d < D; ++d)
  //     r.elts[d] = x.elts[d] + a;
  //   return r;
  // }
  // friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator-(const vect &x,
  //                                                              T a) {
  //   vect r;
  //   for (int d = 0; d < D; ++d)
  //     r.elts[d] = x.elts[d] - a;
  //   return r;
  // }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator*(const T &a,
                                                               const vect &x) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = a * x.elts[d];
    return r;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator*(const vect &x,
                                                               const T &a) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] * a;
    return r;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator/(const vect &x,
                                                               const T &a) {
    vect r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] / a;
    return r;
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator+=(const vect &x) {
    return *this = *this + x;
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator-=(const vect &x) {
    return *this = *this - x;
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator*=(const vect &x) {
    return *this = *this * x;
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator/=(const vect &x) {
    return *this = *this / x;
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator+=(const T &a) {
    return *this = *this + a;
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator-=(const T &a) {
    return *this = *this - a;
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator*=(const T &a) {
    return *this = *this * a;
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect operator/=(const T &a) {
    return *this = *this / a;
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator!(const vect &x) {
    vect<bool, D> r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = !x.elts[d];
    return r;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator&&(const vect &x, const vect &y) {
    vect<bool, D> r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] && y.elts[d];
    return r;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator||(const vect &x, const vect &y) {
    vect<bool, D> r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] || y.elts[d];
    return r;
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator==(const vect &x, const vect &y) {
    vect<bool, D> r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] == y.elts[d];
    return r;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator!=(const vect &x, const vect &y) {
    return !(x == y);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator<(const vect &x, const vect &y) {
    vect<bool, D> r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = x.elts[d] < y.elts[d];
    return r;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator>(const vect &x, const vect &y) {
    return y < x;
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator<=(const vect &x, const vect &y) {
    return !(x > y);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<bool, D>
  operator>=(const vect &x, const vect &y) {
    return !(x < y);
  }

  template <typename U>
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<U, D>
  ifelse(const Arith::vect<U, D> &x, const Arith::vect<U, D> &y) const {
    vect<U, D> r;
    for (int d = 0; d < D; ++d)
      r.elts[d] = elts[d] ? x.elts[d] : y.elts[d];
    return r;
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE bool all() const {
    bool r = true;
    for (int d = 0; d < D; ++d)
      r &= elts[d];
    return r;
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE bool any() const {
    bool r = false;
    for (int d = 0; d < D; ++d)
      r |= elts[d];
    return r;
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE T maxabs() const {
    T r = zero<T>()();
    for (int d = 0; d < D; ++d)
      r = fmax(r, fabs(elts[d]));
    return r;
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE T prod() const {
    T r = 1;
    for (int d = 0; d < D; ++d)
      r *= elts[d];
    return r;
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE T sum() const {
    T r = zero<T>()();
    for (int d = 0; d < D; ++d)
      r += elts[d];
    return r;
  }

  friend ostream &operator<<(ostream &os, const vect &x) {
    os << "[";
    for (int d = 0; d < D; ++d) {
      if (d > 0)
        os << ",";
      os << x.elts[d];
    }
    os << "]";
    return os;
  }
};

template <typename T, int D> struct zero<vect<T, D> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE vect<T, D> operator()() {
    return vect<T, D>();
  }
};

} // namespace Arith
namespace std {
template <typename T, int D> struct equal_to<Arith::vect<T, D> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE bool
  operator()(const Arith::vect<T, D> &lhs, const Arith::vect<T, D> &rhs) const {
    // This is not yet constexpr in C++17
    // return equal_to<array<T, D> >()(lhs.elts, rhs.elts);
    for (int d = 0; d < D; ++d)
      if (!equal_to<T>()(lhs.elts[d], rhs.elts[d]))
        return false;
    return true;
  }
};

template <typename T, int D> struct less<Arith::vect<T, D> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE bool
  operator()(const Arith::vect<T, D> &lhs, const Arith::vect<T, D> &rhs) const {
    return less<array<T, D> >(lhs.elts, rhs.elts);
  }
};
} // namespace std
namespace Arith {

template <typename T, int D>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE Arith::vect<T, D>
abs(const Arith::vect<T, D> &x) {
  Arith::vect<T, D> r;
  for (int d = 0; d < D; ++d)
    r.elts[d] = abs(x.elts[d]);
  return r;
}

template <typename T, int D>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE Arith::vect<bool, D>
isnan1(const Arith::vect<T, D> &x) {
  using std::isnan1;
  Arith::vect<bool, D> r;
  for (int d = 0; d < D; ++d)
    r.elts[d] = isnan1(x.elts[d]);
  return r;
}

template <typename T, int D>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE Arith::vect<T, D>
max(const Arith::vect<T, D> &x, const Arith::vect<T, D> &y) {
  using std::max;
  Arith::vect<T, D> r;
  for (int d = 0; d < D; ++d)
    r.elts[d] = max(x.elts[d], y.elts[d]);
  return r;
}

template <typename T, int D>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE Arith::vect<T, D>
min(const Arith::vect<T, D> &x, const Arith::vect<T, D> &y) {
  using std::min;
  Arith::vect<T, D> r;
  for (int d = 0; d < D; ++d)
    r.elts[d] = min(x.elts[d], y.elts[d]);
  return r;
}

} // namespace Arith

#endif // #ifndef VECT_HXX
