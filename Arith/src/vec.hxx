#ifndef VEC_HXX
#define VEC_HXX

#include "defs.hxx"
#include "dual.hxx"
#include "simd.hxx"
#include "vect.hxx"

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>

#include <array>
#include <cassert>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <tuple>
#include <vector>

namespace Arith {
using namespace std;

// Co- and contravariance

enum class dnup_t : bool { dn, up };
constexpr dnup_t DN = dnup_t::dn;
constexpr dnup_t UP = dnup_t::up;
constexpr dnup_t operator!(const dnup_t dnup) { return dnup_t(!bool(dnup)); }
inline ostream &operator<<(ostream &os, const dnup_t dnup) {
  return os << (dnup == DN ? "d" : "u");
}

// Symmetries

enum class symm_t : int { full, symm, anti };
constexpr symm_t FULL = symm_t::full;
constexpr symm_t SYMM = symm_t::symm;
constexpr symm_t ANTI = symm_t::anti;
template <symm_t symm, typename T1, typename T2, typename T3>
constexpr auto if_symm(const T1 &if_full, const T2 &if_symm,
                       const T3 &if_anti) {
  if constexpr (symm == FULL)
    return if_full;
  if constexpr (symm == SYMM)
    return if_symm;
  if constexpr (symm == ANTI)
    return if_anti;
}
template <typename T1, typename T2, typename T3>
constexpr auto if_symm(symm_t symm, const T1 &if_full, const T2 &if_symm,
                       const T3 &if_anti) {
  if (symm == FULL)
    return if_full;
  if (symm == SYMM)
    return if_symm;
  if (symm == ANTI)
    return if_anti;
  abort();
}
inline ostream &operator<<(ostream &os, const symm_t symm) {
  return os << if_symm(symm, "", "s", "a");
}

////////////////////////////////////////////////////////////////////////////////

// Vector
template <typename T, int D, dnup_t dnup> struct vec {

  // template <typename, int, dnup_t> friend class vec;

  constexpr static int N = D;
  vect<T, N> elts;

public:
  typedef T value_type;
  typedef int size_type;
  static constexpr int size_value = N;

private:
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < D);
#endif
    return n;
  }

public:
  // initializes all elements to nan
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec()
      : elts(nan<vect<T, N> >()) {}

  constexpr ARITH_INLINE vec(const vec &) = default;
  constexpr ARITH_INLINE vec(vec &&) = default;
  constexpr ARITH_INLINE vec &operator=(const vec &) = default;
  constexpr ARITH_INLINE vec &operator=(vec &&) = default;

  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(vec<U, D, dnup> x)
      : elts(move(x.elts)) {}

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(initializer_list<T> x)
      : elts(x) {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(vect<T, N> elts)
      : elts(move(elts)) {}
  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(array<T, N> x)
      : elts(move(x)) {}
  // explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(const vector<T>
  // &x) : elts(x) {} explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  // vec(vector<T> &&x) : elts(move(x)) {}

  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  operator vect<T, N>() const {
    return elts;
  }

  template <typename F>
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void loop(const F &f) {
    vect<T, D>::loop(f);
  }
  template <typename F, typename = result_of_t<F(int)> >
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec(const F &f)
      : elts(vect<T, D>::make(f)) {}

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec pure(const T &a) {
    return vect<T, N>::pure(a);
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec unit(const int i) {
    return vect<T, N>::unit(i);
  }

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<int, D, dnup>
  iota() {
    return {vect<int, N>::iota()};
  }

  template <typename F, typename... Args,
            typename R =
                remove_cv_t<remove_reference_t<result_of_t<F(T, Args...)> > > >
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<R, D, dnup>
  fmap(const F &f, const vec &x, const vec<Args, D, dnup> &...args) {
    return fmap(f, x.elts, args.elts...);
  }
  template <typename F, typename... Args>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  fmap_(const F &f, const vec &x, const vec<Args, D, dnup> &...args) {
    fmap_(f, x.elts, args.elts...);
  }

  template <
      typename... Args,
      typename R = remove_cv_t<remove_reference_t<result_of_t<T(Args...)> > > >
  ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<R, D, dnup>
  operator()(const Args &...args) const {
    return fmap([&](const T &var) { return var(args...); }, *this);
  }
  template <typename Arg1, typename Arg2, typename U,
            typename R = remove_cv_t<
                remove_reference_t<result_of_t<T(Arg1, Arg2, U)> > > >
  ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<R, D, dnup>
  operator()(const Arg1 &arg1, const Arg2 &arg2,
             const vec<U, D, dnup> &val) const {
    return fmap([&](const T &var, const U &x) { return var(arg1, arg2, x); },
                *this, val);
  }
  // template <typename... Args, typename U, typename T1 = T,
  //           typename = decltype(declval<T1>().store(declval<Args>()...,
  //                                                   declval<U>()))>
  // ARITH_INLINE ARITH_DEVICE ARITH_HOST void store(const Args &...args,
  //                         const vec<U, D, dnup> &val) const {
  //   fmap_([&](const auto &var, const auto &x) { return var.store(args..., x);
  //   },
  //        val);
  // }
  template <typename Arg1, typename Arg2, typename U, typename T1 = T,
            typename = decltype(declval<T1>().store(
                declval<Arg1>(), declval<Arg2>(), declval<U>()))>
  ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  store(const Arg1 &arg1, const Arg2 &arg2, const vec<U, D, dnup> &val) const {
    fmap_([&](const T &var, const U &x) { return var.store(arg1, arg2, x); },
          *this, val);
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST size_type size() const {
    return N;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const T &
  operator()(int i) const {
    return elts[ind(i)];
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T &operator()(int i) {
    return elts[ind(i)];
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D, dnup>
  operator+(const vec<T, D, dnup> &x) {
    return {+x.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D, dnup>
  operator-(const vec<T, D, dnup> &x) {
    return {-x.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D, dnup>
  operator+(const vec<T, D, dnup> &x, const vec<T, D, dnup> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D, dnup>
  operator-(const vec<T, D, dnup> &x, const vec<T, D, dnup> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D, dnup>
  operator*(const T &a, const vec<T, D, dnup> &x) {
    return {a * x.elts};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D, dnup>
  operator*(const vec<T, D, dnup> &x, const T &a) {
    return {x.elts * a};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D, dnup>
  operator/(const vec<T, D, dnup> &x, const T &a) {
    return {x.elts / a};
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec<T, D, dnup>
  operator%(const vec<T, D, dnup> &x, const T &a) {
    return {x.elts % a};
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator+=(const vec &x) {
    return *this = *this + x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator-=(const vec &x) {
    return *this = *this - x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator*=(const T &a) {
    return *this = *this * a;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator/=(const T &a) {
    return *this = *this / a;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST vec operator%=(const T &a) {
    return *this = *this % a;
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator==(const vec<T, D, dnup> &x, const vec<T, D, dnup> &y) {
    return all(x.elts == y.elts);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator!=(const vec<T, D, dnup> &x, const vec<T, D, dnup> &y) {
    return !(x == y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  all(const vec &x) {
    return all(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  any(const vec &x) {
    return any(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  anyisnan(const vec &x) {
    return anyisnan(x.elts);
  }

  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*vec<bool, D, dnup>*/
      isnan(const vec &x) {
    return isnan(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T maxabs(const vec &x) {
    return maxabs(x.elts);
  }

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

template <typename T, int D, dnup_t dnup> struct zero<vec<T, D, dnup> > {
  typedef vec<T, D, dnup> value_type;
  // static constexpr value_type value = vec<T, D, dnup>::pure(zero<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return vec<T, D, dnup>::pure(zero<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return vec<T, D, dnup>::pure(zero<T>());
  }
};

template <typename T, int D, dnup_t dnup> struct nan<vec<T, D, dnup> > {
  typedef vec<T, D, dnup> value_type;
  // static constexpr value_type value = vec<T, D, dnup>::pure(nan<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return vec<T, D, dnup>::pure(nan<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return vec<T, D, dnup>::pure(nan<T>());
  }
};

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D, dnup_t dnup>
constexpr vec<simd<T>, D, dnup> if_else(const simdl<T> &cond,
                                        const vec<simd<T>, D, dnup> &x,
                                        const vec<simd<T>, D, dnup> &y) {
  return fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
              x, y);
}

template <typename T, typename U, int D, dnup_t dnup>
constexpr vec<dual<simd<T>, U>, D, dnup>
if_else(const simdl<T> &cond, const vec<dual<simd<T>, U>, D, dnup> &x,
        const vec<dual<simd<T>, U>, D, dnup> &y) {
  return fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
              x, y);
}

} // namespace Arith

#endif // #ifndef VEC_HXX
