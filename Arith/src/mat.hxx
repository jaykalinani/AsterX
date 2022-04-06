#ifndef MAT_HXX
#define MAT_HXX

#include "defs.hxx"
#include "simd.hxx"
#include "vect.hxx"

#include "vec.hxx" // for dnup_t, symm_t

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <utility>
#include <vector>

namespace Arith {
using namespace std;

// A matrix, i.e. a rank-2 tensor. It can have no symmetry (full
// storage), or be symmetric, or anti-symmetric.
template <typename T, int D, dnup_t dnup1, dnup_t dnup2, symm_t symm>
struct gmat {
  static_assert(symm == symm_t::full || dnup1 == dnup2, "");

  constexpr static int Nfull = D * D;
  constexpr static int Nsymm = D * (D + 1) / 2;
  constexpr static int Nanti = D * (D - 1) / 2;
  constexpr static int N = if_symm(symm, Nfull, Nsymm, Nanti);
  vect<T, N> elts;

public:
  typedef T value_type;
  typedef int size_type;
  static constexpr int size_value = N;

private:
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int
  ind_full(const int i, const int j) {
    return i * D + j;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int
  ind_symm(const int i0, const int j0) {
    using std::max, std::min;
    const int i = min(i0, j0);
    const int j = max(i0, j0);
    return i * (2 * D - 1 - i) / 2 + j;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int
  ind_anti(const int i0, const int j0) {
    using std::max, std::min;
    if (i0 == j0)
      return 0;
    const int i = min(i0, j0);
    const int j = max(i0, j0);
    return i * (2 * D - 3 - i) / 2 + j - 1;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int ind(const int i,
                                                                const int j) {
#ifdef CCTK_DEBUG
    assert(i >= 0 && i < D && j >= 0 && j < D);
#endif
    int n{};
    if constexpr (symm == FULL)
      n = ind_full(i, j);
    if constexpr (symm == SYMM)
      n = ind_symm(i, j);
    if constexpr (symm == ANTI)
      n = ind_anti(i, j);
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < N);
#endif
    return n;
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST int sign(const int i,
                                                                 const int j) {
    if constexpr (symm != ANTI)
      return 1;
    if (i == j)
      return 0;
    return i < j ? 1 : -1;
  }

public:
  // initializes all elements to nan
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat()
      : elts(nan<vect<T, N> >()) {}

  constexpr ARITH_INLINE gmat(const gmat &) = default;
  constexpr ARITH_INLINE gmat(gmat &&) = default;
  constexpr ARITH_INLINE gmat &operator=(const gmat &) = default;
  constexpr ARITH_INLINE gmat &operator=(gmat &&) = default;

  template <typename U>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  gmat(gmat<U, D, dnup1, dnup2, symm> x)
      : elts(move(x.elts)) {}

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat(initializer_list<T> x)
      : elts(x) {}
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat(vect<T, N> elts)
      : elts(move(elts)) {}
  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat(array<T, N> x)
      : elts(move(x)) {}
  // explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat(const
  // vector<T> &x) : elts(x) {} explicit constexpr ARITH_INLINE ARITH_DEVICE
  // ARITH_HOST gmat(vector<T> &&x) : elts(move(x)) {}

  explicit constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST
  operator vect<T, N>() const {
    return elts;
  }

  template <typename F>
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void loop(const F &f) {
    if constexpr (symm == symm_t::full)
      for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
          f(i, j);
    if constexpr (symm == symm_t::symm)
      for (int i = 0; i < D; ++i)
        for (int j = i; j < D; ++j)
          f(i, j);
    if constexpr (symm == symm_t::anti)
      for (int i = 0; i < D; ++i)
        for (int j = i + 1; j < D; ++j)
          f(i, j);
  }
  template <typename F, typename = result_of_t<F(int, int)> >
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat(const F &f)
      : gmat(fmap(f, iota1(), iota2())) {}

  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat pure(const T &a) {
    return {vect<T, N>::pure(a)};
  }
  static constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat unit(const int i,
                                                                  const int j) {
    gmat r = zero<gmat>();
    r(i, j) = one<T>();
    return r;
  }

  static constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<int, D, dnup1, dnup2, symm>
      iota() {
    gmat<int, D, dnup1, dnup2, symm> r;
    gmat<int, D, dnup1, dnup2, symm>::loop([&](int i, int j) {
      r(i, j) = {i, j};
    });
    return r;
  }
  static constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<int, D, dnup1, dnup2, symm>
      iota1() {
    gmat<int, D, dnup1, dnup2, symm> r;
    gmat<int, D, dnup1, dnup2, symm>::loop(
        [&](int i, int j) { r.set(i, j, i); });
    return r;
  }
  static constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<int, D, dnup1, dnup2, symm>
      iota2() {
    gmat<int, D, dnup1, dnup2, symm> r;
    gmat<int, D, dnup1, dnup2, symm>::loop(
        [&](int i, int j) { r.set(i, j, j); });
    return r;
  }

  template <typename F, typename... Args,
            typename R =
                remove_cv_t<remove_reference_t<result_of_t<F(T, Args...)> > > >
  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<R, D, dnup1, dnup2, symm>
      fmap(const F &f, const gmat &x,
           const gmat<Args, D, dnup1, dnup2, symm> &...args) {
    return fmap(f, x.elts, args.elts...);
  }
  template <typename F, typename... Args>
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  fmap_(const F &f, const gmat &x,
        const gmat<Args, D, dnup1, dnup2, symm> &...args) {
    fmap_(f, x.elts, args.elts...);
  }

  template <
      typename... Args,
      typename R = remove_cv_t<remove_reference_t<result_of_t<T(Args...)> > > >
  ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<R, D, dnup1, dnup2, symm>
  operator()(const Args &...args) const {
    return fmap([&](const T &var) { return var(args...); }, *this);
  }
  template <typename Arg1, typename Arg2, typename U, typename T1 = T,
            typename R = remove_cv_t<
                remove_reference_t<result_of_t<T(Arg1, Arg2, U)> > > >
  ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<R, D, dnup1, dnup2, symm>
  operator()(const Arg1 &arg1, const Arg2 &arg2,
             const gmat<U, D, dnup1, dnup2, symm> &val) const {
    return fmap([&](const T &var, const U &x) { return var(arg1, arg2, x); },
                *this, val);
  }
  // template <typename... Args, typename U, typename T1 = T,
  //           typename = decltype(declval<T1>().store(declval<Args>()...,
  //                                                   declval<U>()))>
  // ARITH_INLINE ARITH_DEVICE ARITH_HOST void store(const Args &...args,
  //                         const gmat<U, D, dnup1, dnup2, symm> &val) const {
  //   fmap_([&](const auto &var, const auto &x) { return var.store(args..., x);
  //   },
  //        val);
  // }
  template <typename Arg1, typename Arg2, typename U, typename T1 = T,
            typename = decltype(declval<T1>().store(
                declval<Arg1>(), declval<Arg2>(), declval<U>()))>
  ARITH_INLINE ARITH_DEVICE ARITH_HOST void
  store(const Arg1 &arg1, const Arg2 &arg2,
        const gmat<U, D, dnup1, dnup2, symm> &val) const {
    fmap_([&](const T &var, const U &x) { return var.store(arg1, arg2, x); },
          *this, val);
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST size_type size() const {
    return N;
  }
  template <symm_t symm1 = symm, enable_if_t<symm1 != ANTI> * = nullptr>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST const T &
  operator()(int i, int j) const {
    return elts[ind(i, j)];
  }
  template <symm_t symm1 = symm, enable_if_t<symm1 == ANTI> * = nullptr>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T operator()(int i,
                                                              int j) const {
    return sign(i, j) * elts[ind(i, j)];
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST void set(int i, int j,
                                                          const T &x) {
#ifdef CCTK_DEBUG
    assert(sign(i, j) == 1);
#endif
    elts[ind(i, j)] = x;
  }
  template <symm_t symm1 = symm, enable_if_t<symm1 != ANTI> * = nullptr>
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T &operator()(int i, int j) {
    return elts[ind(i, j)];
  }

  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<T, D, dnup1, dnup2, symm>
      operator+(const gmat<T, D, dnup1, dnup2, symm> &x) {
    return {+x.elts};
  }
  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<T, D, dnup1, dnup2, symm>
      operator-(const gmat<T, D, dnup1, dnup2, symm> &x) {
    return {-x.elts};
  }
  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<T, D, dnup1, dnup2, symm>
      operator+(const gmat<T, D, dnup1, dnup2, symm> &x,
                const gmat<T, D, dnup1, dnup2, symm> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<T, D, dnup1, dnup2, symm>
      operator-(const gmat<T, D, dnup1, dnup2, symm> &x,
                const gmat<T, D, dnup1, dnup2, symm> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<T, D, dnup1, dnup2, symm>
      operator*(const T &a, const gmat<T, D, dnup1, dnup2, symm> &x) {
    return {a * x.elts};
  }
  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<T, D, dnup1, dnup2, symm>
      operator*(const gmat<T, D, dnup1, dnup2, symm> &x, const T &a) {
    return {x.elts * a};
  }
  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<T, D, dnup1, dnup2, symm>
      operator/(const gmat<T, D, dnup1, dnup2, symm> &x, const T &a) {
    return {x.elts / a};
  }
  friend constexpr
      ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat<T, D, dnup1, dnup2, symm>
      operator%(const gmat<T, D, dnup1, dnup2, symm> &x, const T &a) {
    return {x.elts % a};
  }

  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat
  operator+=(const gmat &x) {
    return *this = *this + x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat
  operator-=(const gmat &x) {
    return *this = *this - x;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat operator*=(const T &a) {
    return *this = *this * a;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat operator/=(const T &a) {
    return *this = *this / a;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST gmat operator%=(const T &a) {
    return *this = *this % a;
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator==(const gmat<T, D, dnup1, dnup2, symm> &x,
             const gmat<T, D, dnup1, dnup2, symm> &y) {
    return all(x.elts == y.elts);
  }
  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  operator!=(const gmat<T, D, dnup1, dnup2, symm> &x,
             const gmat<T, D, dnup1, dnup2, symm> &y) {
    return !(x == y);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  all(const gmat &x) {
    return all(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST auto /*bool*/
  any(const gmat &x) {
    return any(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST bool
  anyisnan(const gmat &x) {
    return anyisnan(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE
      ARITH_HOST auto /*gmat<bool, D, dnup1, dnup2, symm>*/
      isnan(const gmat &x) {
    return isnan(x.elts);
  }

  friend constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST T
  maxabs(const gmat &x) {
    return maxabs(x.elts);
  }

  friend ostream &operator<<(ostream &os,
                             const gmat<T, D, dnup1, dnup2, symm> &x) {
    os << "(" << dnup1 << dnup2 << symm << ")[";
    for (int j = 0; j < D; ++j) {
      if (j > 0)
        os << ",";
      os << "[";
      for (int i = 0; i < D; ++i) {
        if (i > 0)
          os << ",";
        os << x(i, j);
      }
      os << "]";
    }
    os << "]";
    return os;
  }
};

template <typename T, int D, dnup_t dnup1, dnup_t dnup2, symm_t symm>
struct zero<gmat<T, D, dnup1, dnup2, symm> > {
  typedef gmat<T, D, dnup1, dnup2, symm> value_type;
  static constexpr value_type value =
      gmat<T, D, dnup1, dnup2, symm>::pure(zero<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return value;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return value;
  }
};

// Full and symmetric gmatrices have a unit, antisymmetric gmatrices don't
template <typename T, int D, dnup_t dnup1, dnup_t dnup2>
struct one<gmat<T, D, dnup1, dnup2, symm_t::full> > {
  typedef gmat<T, D, dnup1, dnup2, symm_t::full> value_type;
  static constexpr value_type value = gmat<T, D, dnup1, dnup2, symm_t::full>(
      [](int i, int j) { return i == j ? one<T>()() : zero<T>()(); });
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return value;
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return value;
  }
};
template <typename T, int D, dnup_t dnup1, dnup_t dnup2>
struct one<gmat<T, D, dnup1, dnup2, symm_t::symm> > {
  typedef gmat<T, D, dnup1, dnup2, symm_t::symm> value_type;
  // static constexpr value_type value = gmat<T, D, dnup1, dnup2, symm_t::symm>(
  //     [](int i, int j) { return i == j; });
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return gmat<T, D, dnup1, dnup2, symm_t::symm>(
        [](int i, int j) { return i == j ? one<T>()() : zero<T>()(); });
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return gmat<T, D, dnup1, dnup2, symm_t::symm>(
        [](int i, int j) { return i == j ? one<T>()() : zero<T>()(); });
  }
};

template <typename T, int D, dnup_t dnup1, dnup_t dnup2, symm_t symm>
struct nan<gmat<T, D, dnup1, dnup2, symm> > {
  typedef gmat<T, D, dnup1, dnup2, symm> value_type;
  // static constexpr value_type value =
  //     gmat<T, D, dnup1, dnup2, symm>::pure(nan<T>::value);
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST operator value_type() const {
    return gmat<T, D, dnup1, dnup2, symm>::pure(nan<T>());
  }
  constexpr ARITH_INLINE ARITH_DEVICE ARITH_HOST value_type operator()() const {
    return gmat<T, D, dnup1, dnup2, symm>::pure(nan<T>());
  }
};

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D, dnup_t dnup1, dnup_t dnup2>
using mat = gmat<T, D, dnup1, dnup2, symm_t::full>;
template <typename T, int D, dnup_t dnup1, dnup_t dnup2>
using smat = gmat<T, D, dnup1, dnup2, symm_t::symm>;
template <typename T, int D, dnup_t dnup1, dnup_t dnup2>
using amat = gmat<T, D, dnup1, dnup2, symm_t::anti>;

////////////////////////////////////////////////////////////////////////////////

template <typename T, int D, dnup_t dnup1, dnup_t dnup2, symm_t symm>
constexpr gmat<simd<T>, D, dnup1, dnup2, symm>
if_else(const simdl<T> &cond, const gmat<simd<T>, D, dnup1, dnup2, symm> &x,
        const gmat<simd<T>, D, dnup1, dnup2, symm> &y) {
  return fmap([&](const auto &x, const auto &y) { return if_else(cond, x, y); },
              x, y);
}

} // namespace Arith

#endif // #ifndef MAT_HXX
