#ifndef MAT_HXX
#define MAT_HXX

#include "vect.hxx"

// Only for dnup_t
#include "vec.hxx"

#include <cctk.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <functional>
#include <initializer_list>
#include <iostream>
#include <tuple>
#include <utility>
#include <vector>

#ifdef CCTK_DEBUG
#define ARITH_INLINE
#else
#define ARITH_INLINE CCTK_ATTRIBUTE_ALWAYS_INLINE
#endif

namespace Arith {

// Symmetric 3-matrix
template <typename T, int D, dnup_t dnup1, dnup_t dnup2> class mat {
  static_assert(dnup1 == dnup2, "");

  template <typename, int, dnup_t, dnup_t> friend class mat;

  constexpr static int N = D * (D - 1);
  vect<T, N> elts;

  static constexpr ARITH_INLINE int symind(const int i, const int j) {
#ifdef CCTK_DEBUG
    assert(i >= 0 && i <= j && j < 3);
#endif
    const int n = 3 * i - i * (i + 1) / 2 + j;
    // i j n
    // 0 0 0
    // 0 1 1
    // 0 2 2
    // 1 1 3
    // 1 2 4
    // 2 2 5
    // const int n = 2 * i + j - (unsigned)i / 2;
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < N);
#endif
    return n;
  }
  static constexpr ARITH_INLINE int ind(const int i, const int j) {
    using std::max, std::min;
    return symind(min(i, j), max(i, j));
  }

  static_assert(symind(0, 0) == 0, "");
  static_assert(symind(0, 1) == 1, "");
  static_assert(symind(0, 2) == 2, "");
  static_assert(symind(1, 1) == 3, "");
  static_assert(symind(1, 2) == 4, "");
  static_assert(symind(2, 2) == 5, "");

  static_assert(ind(1, 0) == ind(0, 1), "");
  static_assert(ind(2, 0) == ind(0, 2), "");
  static_assert(ind(2, 1) == ind(1, 2), "");

public:
  explicit constexpr ARITH_INLINE mat() : elts() {}

  constexpr ARITH_INLINE mat(const mat &) = default;
  constexpr ARITH_INLINE mat(mat &&) = default;
  constexpr ARITH_INLINE mat &operator=(const mat &) = default;
  constexpr ARITH_INLINE mat &operator=(mat &&) = default;

  template <typename U>
  constexpr ARITH_INLINE mat(const mat<U, D, dnup1, dnup2> &x) : elts(x.elts) {}
  template <typename U>
  constexpr ARITH_INLINE mat(mat<U, D, dnup1, dnup2> &&x)
      : elts(move(x.elts)) {}

  constexpr ARITH_INLINE mat(const vect<T, N> &elts) : elts(elts) {}
  constexpr ARITH_INLINE mat(vect<T, N> &&elts) : elts(move(elts)) {}

  constexpr ARITH_INLINE mat(initializer_list<T> A) : elts(A) {}
  constexpr ARITH_INLINE mat(const array<T, 6> &A) : elts(A) {}
  constexpr ARITH_INLINE mat(array<T, 6> &&A) : elts(move(A)) {}
  constexpr ARITH_INLINE mat(const vector<T> &A) : elts(A) {}
  constexpr ARITH_INLINE mat(vector<T> &&A) : elts(move(A)) {}

  template <typename F, typename = result_of_t<F(int, int)> >
  constexpr ARITH_INLINE mat(F f) : mat(iota1().map(f, iota2())) {
    // #ifdef CCTK_DEBUG
    //     // Check symmetry
    //     const auto is_symmetric{[](const T &fgood, const T &fother) {
    //       return norm1<T>()(fother - fgood) <=
    //              1.0e-12 * (1 + norm1<T>()(fgood) + norm1<T>()(fother));
    //     }};
    //     for (int i = 0; i < D; ++i) {
    //       for (j = i + 1; j < D; ++j) {
    //         const T fij = f(i, j);
    //         const T fji = f(j, i);
    //         if (!is_symmetric(fij, fji)) {
    //           ostringstream buf;
    //           buf << "f(0,1)=" << f(0, 1) << "\n"
    //               << "f(1,0)=" << f10 << "\n"
    //               << "f(0,2)=" << f(0, 2) << "\n"
    //               << "f(2,0)=" << f20 << "\n"
    //               << "f(1,2)=" << f(1, 2) << "\n"
    //               << "f(2,1)=" << f21 << "\n";
    //           CCTK_VERROR("symmetric matrix is not symmetric:\n%s",
    //                       buf.str().c_str());
    //         }
    //         assert(is_symmetric(fij));
    //       }
    //     }
    // #endif
  }

  static constexpr ARITH_INLINE mat unit(int i, int j) {
    mat r;
    r(i, j) = 1;
    return r;
  }

  static constexpr ARITH_INLINE mat<array<int, 2>, D, dnup1, dnup2> iota() {
    mat<array<int, 2>, D, dnup1, dnup2> r;
    for (int i = 0; i < D; ++i)
      for (int j = i; j < D; ++j)
        r(i, j) = {i, j};
    return r;
  }
  static constexpr ARITH_INLINE mat<int, D, dnup1, dnup2> iota1() {
    mat<int, D, dnup1, dnup2> r;
    for (int i = 0; i < D; ++i)
      for (int j = i; j < D; ++j)
        r(i, j) = i;
    return r;
  }
  static constexpr ARITH_INLINE mat<int, D, dnup1, dnup2> iota2() {
    mat<int, D, dnup1, dnup2> r;
    for (int i = 0; i < D; ++i)
      for (int j = i; j < D; ++j)
        r(i, j) = j;
    return r;
  }

  template <typename F,
            typename R = remove_cv_t<remove_reference_t<result_of_t<F(T)> > > >
  constexpr ARITH_INLINE mat<R, D, dnup1, dnup2> map(F f) const {
    mat<R, D, dnup1, dnup2> r;
    for (int i = 0; i < N; ++i)
      r.elts[i] = f(elts[i]);
    return r;
  }
  template <
      typename F, typename U,
      typename R = remove_cv_t<remove_reference_t<result_of_t<F(T, U)> > > >
  constexpr ARITH_INLINE mat<R, D, dnup1, dnup2>
  map(F f, const mat<U, D, dnup1, dnup2> &x) const {
    mat<R, D, dnup1, dnup2> r;
    for (int i = 0; i < N; ++i)
      r.elts[i] = f(elts[i], x.elts[i]);
    return r;
  }

  constexpr ARITH_INLINE const T &operator()(int i, int j) const {
    return elts[ind(i, j)];
  }
  constexpr ARITH_INLINE T &operator()(int i, int j) { return elts[ind(i, j)]; }

  // template <typename U = T>
  // ARITH_INLINE
  //     mat3<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >,
  //          dnup1, dnup2>
  //     operator()(const vect<int, 3> &I) const {
  //   return {elts[0](I), elts[1](I), elts[2](I),
  //           elts[3](I), elts[4](I), elts[5](I)};
  // }

  friend constexpr ARITH_INLINE mat<T, D, dnup1, dnup2>
  operator+(const mat<T, D, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend constexpr ARITH_INLINE mat<T, D, dnup1, dnup2>
  operator-(const mat<T, D, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend constexpr ARITH_INLINE mat<T, D, dnup1, dnup2>
  operator+(const mat<T, D, dnup1, dnup2> &x,
            const mat<T, D, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr ARITH_INLINE mat<T, D, dnup1, dnup2>
  operator-(const mat<T, D, dnup1, dnup2> &x,
            const mat<T, D, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr ARITH_INLINE mat<T, D, dnup1, dnup2>
  operator*(const T &a, const mat<T, D, dnup1, dnup2> &x) {
    return {a * x.elts};
  }
  friend constexpr ARITH_INLINE mat<T, D, dnup1, dnup2>
  operator*(const mat<T, D, dnup1, dnup2> &x, const T &a) {
    return {x.elts * a};
  }

  constexpr ARITH_INLINE mat operator+=(const mat &x) {
    return *this = *this + x;
  }
  constexpr ARITH_INLINE mat operator-=(const mat &x) {
    return *this = *this - x;
  }
  constexpr ARITH_INLINE mat operator*=(const T &a) {
    return *this = *this * a;
  }
  constexpr ARITH_INLINE mat operator/=(const T &a) {
    return *this = *this / a;
  }

  friend constexpr ARITH_INLINE bool
  operator==(const mat<T, D, dnup1, dnup2> &x,
             const mat<T, D, dnup1, dnup2> &y) {
    return equal_to<vect<T, N> >()(x.elts, y.elts);
  }
  friend constexpr ARITH_INLINE bool
  operator!=(const mat<T, D, dnup1, dnup2> &x,
             const mat<T, D, dnup1, dnup2> &y) {
    return !(x == y);
  }

  constexpr ARITH_INLINE T maxabs() const { return elts.maxabs(); }

  // constexpr ARITH_INLINE T det() const {
  //   const auto &A = *this;
  //   return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) -
  //          A(1, 0) * (A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) +
  //          A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
  // }

  // constexpr ARITH_INLINE mat3<T, !dnup1, !dnup2> inv(const T detA) const {
  //   const auto &A = *this;
  //   const T detA1 = 1 / detA;
  //   return mat3<T, !dnup1, !dnup2>{
  //       detA1 * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)),
  //       detA1 * (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)),
  //       detA1 * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)),
  //       detA1 * (A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)),
  //       detA1 * (A(2, 0) * A(0, 1) - A(2, 1) * A(0, 0)),
  //       detA1 * (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0))};
  // }

  // constexpr ARITH_INLINE T trace(const mat3<T, !dnup1, !dnup2> &gu) const {
  //   const auto &A = *this;
  //   return sum2([&](int x, int y) ARITH_INLINE { return gu(x, y) * A(x, y);
  //   });
  // }

  // constexpr ARITH_INLINE mat3 trace_free(
  //     const mat<T, D,dnup1, dnup2> &g, const mat3<T, !dnup1, !dnup2> &gu)
  //     const {
  //   const auto &A = *this;
  //   const T trA = A.trace(gu);
  //   return mat3([&](int a, int b)
  //                   ARITH_INLINE { return A(a, b) - trA / 3 * g(a, b); });
  // }

  friend ostream &operator<<(ostream &os, const mat<T, D, dnup1, dnup2> &A) {
    os << "(" << dnup1 << dnup2 << ")[";
    for (int j = 0; j < D; ++j) {
      if (j > 0)
        os << ",";
      os << "[";
      for (int i = 0; i < D; ++i) {
        if (i > 0)
          os << ",";
        os << A(i, j);
      }
      os << "]";
    }
    os << "]";
    return os;
  }
}; // namespace Arith

} // namespace Arith

#endif // #ifndef MAT_HXX
