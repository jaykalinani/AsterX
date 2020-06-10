#ifndef VEC_HXX
#define VEC_HXX

#include "vect.hxx"

#include <cctk.h>

#include <array>
#include <cassert>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <vector>

#ifdef CCTK_DEBUG
#define ARITH_INLINE
#else
#define ARITH_INLINE CCTK_ATTRIBUTE_ALWAYS_INLINE
#endif

namespace Arith {
using namespace std;

enum class dnup_t : bool { dn, up };
constexpr dnup_t DN = dnup_t::dn;
constexpr dnup_t UP = dnup_t::up;
constexpr dnup_t operator!(const dnup_t dnup) { return dnup_t(!bool(dnup)); }
inline ostream &operator<<(ostream &os, const dnup_t dnup) {
  return os << (dnup == DN ? "d" : "u");
}

// vector
template <typename T, int D, dnup_t dnup> class vec {

  template <typename, int, dnup_t> friend class vec;

  vect<T, D> elts;

  static constexpr ARITH_INLINE int ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < D);
#endif
    return n;
  }

public:
  // initializes all elements to zero
  constexpr ARITH_INLINE vec() : elts() {}

  constexpr ARITH_INLINE vec(const vec &) = default;
  constexpr ARITH_INLINE vec(vec &&) = default;
  constexpr ARITH_INLINE vec &operator=(const vec &) = default;
  constexpr ARITH_INLINE vec &operator=(vec &&) = default;

  template <typename U>
  constexpr ARITH_INLINE vec(const vec<U, D, dnup> &x) : elts(x.elts) {}
  template <typename U>
  constexpr ARITH_INLINE vec(vec<U, D, dnup> &&x) : elts(move(x.elts)) {}

  constexpr ARITH_INLINE vec(const vect<T, D> &elts) : elts(elts) {}
  constexpr ARITH_INLINE vec(vect<T, D> &&elts) : elts(move(elts)) {}

  constexpr ARITH_INLINE vec(initializer_list<T> v) : elts(v) {}
  constexpr ARITH_INLINE vec(const array<T, D> &v) : elts(v) {}
  constexpr ARITH_INLINE vec(array<T, D> &&v) : elts(move(v)) {}
  constexpr ARITH_INLINE vec(const vector<T> &v) : elts(v) {}
  constexpr ARITH_INLINE vec(vector<T> &&v) : elts(move(v)) {}

  template <typename F, typename = result_of_t<F(int)> >
  constexpr ARITH_INLINE vec(F f) : vec(iota().map(f)) {}

  static constexpr ARITH_INLINE vec unit(int i) {
    vec r;
    r(i) = 1;
    return r;
  }

  static constexpr ARITH_INLINE vec<int, D, dnup> iota() {
    vec<int, D, dnup> r;
    for (int i = 0; i < D; ++i)
      r(i) = i;
    return r;
  }

  template <typename F,
            typename R = remove_cv_t<remove_reference_t<result_of_t<F(T)> > > >
  constexpr ARITH_INLINE vec<R, D, dnup> map(F f) const {
    return vec<R, D, dnup>(elts.map(f));
  }
  template <
      typename F, typename U,
      typename R = remove_cv_t<remove_reference_t<result_of_t<F(T, U)> > > >
  constexpr ARITH_INLINE vec<R, D, dnup> map(F f,
                                             const vec<U, D, dnup> &x) const {
    return vec<R, D, dnup>(elts.map(f, x.elt));
  }

  constexpr ARITH_INLINE const T &operator()(int i) const {
    return elts[ind(i)];
  }
  constexpr ARITH_INLINE T &operator()(int i) { return elts[ind(i)]; }

  // template <typename U = T>
  // ARITH_INLINE
  //     vec3<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >,
  //          dnup>
  //     operator()(const vect<int, 3> &I) const {
  //   return {elts[0](I), elts[1](I), elts[2](I)};
  // }

  friend constexpr ARITH_INLINE vec<T, D, dnup>
  operator+(const vec<T, D, dnup> &x) {
    return {+x.elts};
  }
  friend constexpr ARITH_INLINE vec<T, D, dnup>
  operator-(const vec<T, D, dnup> &x) {
    return {-x.elts};
  }
  friend constexpr ARITH_INLINE vec<T, D, dnup>
  operator+(const vec<T, D, dnup> &x, const vec<T, D, dnup> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr ARITH_INLINE vec<T, D, dnup>
  operator-(const vec<T, D, dnup> &x, const vec<T, D, dnup> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr ARITH_INLINE vec<T, D, dnup>
  operator*(const T &a, const vec<T, D, dnup> &x) {
    return {a * x.elts};
  }
  friend constexpr ARITH_INLINE vec<T, D, dnup>
  operator*(const vec<T, D, dnup> &x, const T &a) {
    return {x.elts * a};
  }

  constexpr ARITH_INLINE vec operator+=(const vec &x) {
    return *this = *this + x;
  }
  constexpr ARITH_INLINE vec operator-=(const vec &x) {
    return *this = *this - x;
  }
  constexpr ARITH_INLINE vec operator*=(const T &a) {
    return *this = *this * a;
  }
  constexpr ARITH_INLINE vec operator/=(const T &a) {
    return *this = *this / a;
  }

  friend constexpr ARITH_INLINE bool operator==(const vec<T, D, dnup> &x,
                                                const vec<T, D, dnup> &y) {
    return equal_to<vect<T, D> >()(x.elts, y.elts);
  }
  friend constexpr ARITH_INLINE bool operator!=(const vec<T, D, dnup> &x,
                                                const vec<T, D, dnup> &y) {
    return !(x == y);
  }

  constexpr ARITH_INLINE T maxabs() const { return elts.maxabs(); }

  friend ostream &operator<<(ostream &os, const vec<T, D, dnup> &v) {
    os << "(" << dnup << ")[";
    for (int i = 0; i < D; ++i) {
      if (i > 0)
        os << ",";
      os << v(i);
    }
    os << "]";
    return os;
  }
};

} // namespace Arith

#undef ARITH_INLINE

#endif // #ifndef VEC_HXX
