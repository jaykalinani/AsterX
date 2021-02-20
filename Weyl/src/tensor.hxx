#ifndef TENSOR_HXX
#define TENSOR_HXX

#include <dual.hxx>
#include <loop.hxx>
#include <vect.hxx>

#include <cctk.h>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <sstream>
#include <type_traits>

namespace Weyl {
using namespace Arith;
using namespace Loop;
using namespace std;

template <typename T> constexpr T pow2(const T x) { return x * x; }
template <typename T> constexpr T pow3(const T x) {
  const T x2 = x * x;
  return x2 * x;
}
template <typename T> constexpr T pow4(const T x) {
  const T x2 = x * x;
  return x2 * x2;
}
template <typename T> constexpr T pow5(const T x) {
  const T x2 = x * x;
  const T x4 = x2 * x2;
  return x4 * x;
}
template <typename T> constexpr T pow6(const T x) {
  const T x2 = x * x;
  const T x4 = x2 * x2;
  return x4 * x2;
}

namespace detail {
template <typename T> constexpr T pown(const T x, int n) {
  T r{1};
  T y{x};
  while (n) {
    if (n & 1)
      r *= y;
    y *= y;
    n >>= 1;
  }
  return r;
}
} // namespace detail

template <typename T> constexpr T pown(const T x, const int n) {
  return n >= 0 ? detail::pown(x, n) : 1 / detail::pown(x, -n);
}

constexpr int factorial(int n) {
  int r{1};
  while (n > 1) {
    r *= n;
    --n;
  }
  return r;
}

////////////////////////////////////////////////////////////////////////////////

template <typename F,
          typename R = remove_cv_t<remove_reference_t<result_of_t<F(int)> > > >
constexpr R sum41(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 4; ++x)
    s += f(x);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int)> > > >
constexpr R sum42(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 4; ++y)
      s += f(x, y);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int, int)> > > >
constexpr R sum43(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 4; ++y)
      for (int z = 0; z < 4; ++z)
        s += f(x, y, z);
  return s;
}

template <typename F, typename R = remove_cv_t<remove_reference_t<
                          result_of_t<F(int, int, int, int)> > > >
constexpr R sum44(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 4; ++y)
      for (int z = 0; z < 4; ++z)
        for (int w = 0; w < 4; ++w)
          s += f(x, y, z, w);
  return s;
}

////////////////////////////////////////////////////////////////////////////////

template <typename F, typename T> constexpr T fold(const F &f, const T &x) {
  return x;
}
template <typename F, typename T, typename... Ts>
constexpr T fold(const F &f, const T &x0, const T &x1, const Ts &...xs) {
  return fold(f, fold(f, x0, x1), xs...);
}

template <typename T> constexpr T add() { return T(0); }
// template <typename T>
//  constexpr  T add(const T &x) {
//   return x;
// }
template <typename T, typename... Ts>
constexpr T add(const T &x, const Ts &...xs) {
  return x + add(xs...);
}

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct nan {
  constexpr T operator()() const { return NAN; }
};
template <typename T, int D> struct nan<vect<T, D> > {
  constexpr vect<T, D> operator()() const {
    return vect<T, D>::pure(nan<T>()());
  }
};

template <typename T> struct norm1 {
  typedef T result_type;
  constexpr result_type operator()(const T &x) const { return abs(x); }
};
template <typename T, int D> struct norm1<vect<T, D> > {
  typedef typename norm1<T>::result_type result_type;
  constexpr result_type operator()(const vect<T, D> &xs) const {
    typedef typename norm1<T>::result_type R;
    R r = zero<R>()();
    for (int d = 0; d < D; ++d)
      r = max(r, norm1<T>()(xs[d]));
    return r;
  }
};

////////////////////////////////////////////////////////////////////////////////

enum class dnup_t : bool { dn, up };
constexpr dnup_t DN = dnup_t::dn;
constexpr dnup_t UP = dnup_t::up;
constexpr dnup_t operator!(const dnup_t dnup) { return dnup_t(!bool(dnup)); }

// 3-vector
template <typename T, dnup_t dnup> class vec3 {
  vect<T, 3> elts;

  static constexpr int ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 3);
#endif
    return n;
  }

public:
  explicit constexpr vec3() : elts{nan<vect<T, 3> >()()} {}

  constexpr vec3(const vect<T, 3> &elts) : elts(elts) {}
  constexpr vec3(vect<T, 3> &&elts) : elts(move(elts)) {}

  explicit constexpr vec3(T vx, T vy, T vz)
      : elts(make_tuple(move(vx), move(vy), move(vz))) {}

  constexpr vec3(initializer_list<T> v) : elts(v) {}
  // constexpr vec3(const vector<T> &v) : elts(v) {}
  // constexpr vec3(vector<T> &&v) : elts(move(v)) {}

  vec3(const GF3D2<add_const_t<T> > &gf_vx_,
       const GF3D2<add_const_t<T> > &gf_vy_,
       const GF3D2<add_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec3{gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  vec3(const GF3D2<remove_const_t<T> > &gf_vx_,
       const GF3D2<remove_const_t<T> > &gf_vy_,
       const GF3D2<remove_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec3{gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  template <typename F, typename = result_of_t<F(int)> >
  constexpr vec3(const F &f) : elts{f(0), f(1), f(2)} {}

  void store(const GF3D2<T> &gf_vx_, const GF3D2<T> &gf_vy_,
             const GF3D2<T> &gf_vz_, const vect<int, 3> &I) const {
    const auto &v = *this;
#ifdef CCTK_DEBUG
    if (!((CCTK_isfinite(v(0))) && (CCTK_isfinite(v(1))) &&
          (CCTK_isfinite(v(2))))) {
      ostringstream buf;
      buf << "v=" << v;
      CCTK_VERROR("nan found: %s", buf.str().c_str());
    }
    assert(CCTK_isfinite(v(0)));
    assert(CCTK_isfinite(v(1)));
    assert(CCTK_isfinite(v(2)));
#endif
    gf_vx_(I) = v(0);
    gf_vy_(I) = v(1);
    gf_vz_(I) = v(2);
  }

  const T &operator()(int i) const { return elts[ind(i)]; }
  T &operator()(int i) { return elts[ind(i)]; }

  template <typename U = T>
  vec3<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >, dnup>
  operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  vec3<remove_cv_t<
           remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
       dnup>
  operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I)};
  }

  friend constexpr vec3<T, dnup> operator+(const vec3<T, dnup> &x) {
    return {+x.elts};
  }
  friend constexpr vec3<T, dnup> operator-(const vec3<T, dnup> &x) {
    return {-x.elts};
  }
  friend constexpr vec3<T, dnup> operator+(const vec3<T, dnup> &x,
                                           const vec3<T, dnup> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr vec3<T, dnup> operator-(const vec3<T, dnup> &x,
                                           const vec3<T, dnup> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr vec3<T, dnup> operator*(const T &a, const vec3<T, dnup> &x) {
    return {a * x.elts};
  }

  friend constexpr bool operator==(const vec3<T, dnup> &x,
                                   const vec3<T, dnup> &y) {
    return equal_to<vect<T, 3> >()(x.elts, y.elts);
  }
  friend constexpr bool operator!=(const vec3<T, dnup> &x,
                                   const vec3<T, dnup> &y) {
    return !(x == y);
  }

  constexpr T maxabs() const { return elts.maxabs(); }

  friend struct norm1<vec3>;

  friend ostream &operator<<(ostream &os, const vec3<T, dnup> &v) {
    return os << "[" << v(0) << "," << v(1) << "," << v(2) << "]";
  }
};

template <typename T, dnup_t dnup> struct nan<vec3<T, dnup> > {
  constexpr vec3<T, dnup> operator()() const { return vec3<T, dnup>(); }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup> struct zero<Weyl::vec3<T, dnup> > {
  constexpr Weyl::vec3<T, dnup> operator()() const {
    return Weyl::vec3<T, dnup>(zero<vect<T, 3> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup> struct norm1<vec3<T, dnup> > {
  typedef typename norm1<vect<T, 3> >::result_type result_type;
  constexpr result_type operator()(const vec3<T, dnup> &x) const {
    return norm1<vect<T, 3> >()(x.elts);
  }
};

////////////////////////////////////////////////////////////////////////////////

// Symmetric 3-matrix
template <typename T, dnup_t dnup1, dnup_t dnup2> class mat3 {
  static_assert(dnup1 == dnup2, "");

  vect<T, 6> elts;

  static constexpr int symind(const int i, const int j) {
#ifdef CCTK_DEBUG
    assert(i >= 0 && i <= j && j < 3);
#endif
    const int n = i * (5 - i) / 2 + j;
    // i j n
    // 0 0 0
    // 0 1 1
    // 0 2 2
    // 1 1 3
    // 1 2 4
    // 2 2 5
    // const int n = 2 * i + j - (unsigned)i / 2;
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 6);
#endif
    return n;
  }
  static constexpr int ind(const int i, const int j) {
    return symind(min(i, j), max(i, j));
  }

  static_assert(symind(0, 0) == 0, "");
  static_assert(symind(0, 1) == 1, "");
  static_assert(symind(0, 2) == 2, "");
  static_assert(symind(1, 1) == 3, "");
  static_assert(symind(1, 2) == 4, "");
  static_assert(symind(2, 2) == 5, "");

  // nvcc doesn't handle these constexpr expressions
#ifndef __CUDACC__
  static_assert(ind(1, 0) == ind(0, 1), "");
  static_assert(ind(2, 0) == ind(0, 2), "");
  static_assert(ind(2, 1) == ind(1, 2), "");
#endif

public:
  explicit constexpr mat3() : elts{nan<vect<T, 6> >()()} {}

  constexpr mat3(const vect<T, 6> &elts) : elts(elts) {}
  constexpr mat3(vect<T, 6> &&elts) : elts(move(elts)) {}

  explicit constexpr mat3(T Axx, T Axy, T Axz, T Ayy, T Ayz, T Azz)
      : elts(make_tuple(move(Axx), move(Axy), move(Axz), move(Ayy), move(Ayz),
                        move(Azz))) {}

  constexpr mat3(initializer_list<T> A) : elts(A) {}
  // constexpr mat3(const vector<T> &A) : elts(A) {}
  // constexpr mat3(vector<T> &&A) : elts(move(A)) {}

  mat3(const GF3D2<add_const_t<T> > &gf_Axx_,
       const GF3D2<add_const_t<T> > &gf_Axy_,
       const GF3D2<add_const_t<T> > &gf_Axz_,
       const GF3D2<add_const_t<T> > &gf_Ayy_,
       const GF3D2<add_const_t<T> > &gf_Ayz_,
       const GF3D2<add_const_t<T> > &gf_Azz_, const vect<int, 3> &I)
      : mat3{gf_Axx_(I), gf_Axy_(I), gf_Axz_(I),
             gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  mat3(const GF3D2<remove_const_t<T> > &gf_Axx_,
       const GF3D2<remove_const_t<T> > &gf_Axy_,
       const GF3D2<remove_const_t<T> > &gf_Axz_,
       const GF3D2<remove_const_t<T> > &gf_Ayy_,
       const GF3D2<remove_const_t<T> > &gf_Ayz_,
       const GF3D2<remove_const_t<T> > &gf_Azz_, const vect<int, 3> &I)
      : mat3{gf_Axx_(I), gf_Axy_(I), gf_Axz_(I),
             gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  template <typename F, typename = result_of_t<F(int, int)> >
  constexpr mat3(const F &f)
      : elts{f(0, 0), f(0, 1), f(0, 2), f(1, 1), f(1, 2), f(2, 2)} {
    // #ifdef CCTK_DEBUG
    //     // Check symmetry
    //     const T f10 = f(1, 0);
    //     const T f20 = f(2, 0);
    //     const T f21 = f(1, 2);
    //     const auto scale = norm1<T>()(elts[0]) + norm1<T>()(elts[1]) +
    //                        norm1<T>()(elts[2]) + norm1<T>()(elts[3]) +
    //                        norm1<T>()(elts[4]) + norm1<T>()(elts[5]);
    //     const auto is_approx{[scale](const T &fgood, const T &fother) {
    //       return norm1<T>()(fother - fgood) <= 1.0e-12 * (1 + scale);
    //     }};
    //     if (!(is_approx((*this)(0, 1), f10) && is_approx((*this)(0, 2), f20)
    //     &&
    //           is_approx((*this)(1, 2), f21))) {
    //       ostringstream buf;
    //       buf << "f(0,1)=" << (*this)(0, 1) << "\n"
    //           << "f(1,0)=" << f10 << "\n"
    //           << "f(0,2)=" << (*this)(0, 2) << "\n"
    //           << "f(2,0)=" << f20 << "\n"
    //           << "f(1,2)=" << (*this)(1, 2) << "\n"
    //           << "f(2,1)=" << f21 << "\n";
    //       CCTK_VERROR("symmetric matrix is not symmetric:\n%s",
    //       buf.str().c_str());
    //     }
    //     assert(is_approx(f10, (*this)(0, 1)));
    //     assert(is_approx(f20, (*this)(0, 2)));
    //     assert(is_approx(f21, (*this)(1, 2)));
    // #endif
  }

  void store(const GF3D2<T> &gf_Axx_, const GF3D2<T> &gf_Axy_,
             const GF3D2<T> &gf_Axz_, const GF3D2<T> &gf_Ayy_,
             const GF3D2<T> &gf_Ayz_, const GF3D2<T> &gf_Azz_,
             const vect<int, 3> &I) const {
    const auto &A = *this;
#ifdef CCTK_DEBUG
    assert(CCTK_isfinite(A(0, 0)));
    assert(CCTK_isfinite(A(0, 1)));
    assert(CCTK_isfinite(A(0, 2)));
    assert(CCTK_isfinite(A(1, 1)));
    assert(CCTK_isfinite(A(1, 2)));
    assert(CCTK_isfinite(A(2, 2)));
#endif
    gf_Axx_(I) = A(0, 0);
    gf_Axy_(I) = A(0, 1);
    gf_Axz_(I) = A(0, 2);
    gf_Ayy_(I) = A(1, 1);
    gf_Ayz_(I) = A(1, 2);
    gf_Azz_(I) = A(2, 2);
  }

  const T &operator()(int i, int j) const { return elts[ind(i, j)]; }
  //  T &operator()(int i, int j) { return
  // elts[symind(i, j)]; }
  T &operator()(int i, int j) { return elts[ind(i, j)]; }

  template <typename U = T>
  mat3<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >, dnup1,
       dnup2>
  operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I),
            elts[3](I), elts[4](I), elts[5](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  mat3<remove_cv_t<
           remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
       dnup1, dnup2>
  operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[3](layout, I), elts[4](layout, I), elts[5](layout, I)};
  }

  friend constexpr mat3<T, dnup1, dnup2>
  operator+(const mat3<T, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend constexpr mat3<T, dnup1, dnup2>
  operator-(const mat3<T, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend constexpr mat3<T, dnup1, dnup2>
  operator+(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr mat3<T, dnup1, dnup2>
  operator-(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr mat3<T, dnup1, dnup2>
  operator*(const T &a, const mat3<T, dnup1, dnup2> &x) {
    return {a * x.elts};
  }

  friend constexpr bool operator==(const mat3<T, dnup1, dnup2> &x,
                                   const mat3<T, dnup1, dnup2> &y) {
    return equal_to<vect<T, 6> >()(x.elts, y.elts);
  }
  friend constexpr bool operator!=(const mat3<T, dnup1, dnup2> &x,
                                   const mat3<T, dnup1, dnup2> &y) {
    return !(x == y);
  }

  constexpr T maxabs() const { return elts.maxabs(); }

  friend struct norm1<mat3>;

  constexpr T det() const {
    const auto &A = *this;
    return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) -
           A(1, 0) * (A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) +
           A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
  }

  constexpr mat3<T, !dnup1, !dnup2> inv(const T detA) const {
    const auto &A = *this;
    const T detA1 = 1 / detA;
    return mat3<T, !dnup1, !dnup2>{
        detA1 * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)),
        detA1 * (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)),
        detA1 * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0)),
        detA1 * (A(2, 2) * A(0, 0) - A(2, 0) * A(0, 2)),
        detA1 * (A(2, 0) * A(0, 1) - A(2, 1) * A(0, 0)),
        detA1 * (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0))};
  }

  constexpr T trace(const mat3<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    return sum2([&](int x, int y) { return gu(x, y) * A(x, y); });
  }

  constexpr mat3 trace_free(const mat3<T, dnup1, dnup2> &g,
                            const mat3<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    const T trA = A.trace(gu);
    return mat3([&](int a, int b) { return A(a, b) - trA / 3 * g(a, b); });
  }

  friend ostream &operator<<(ostream &os, const mat3<T, dnup1, dnup2> &A) {
    return os << "[[" << A(0, 0) << "," << A(0, 1) << "," << A(0, 2) << "],["
              << A(1, 0) << "," << A(1, 1) << "," << A(1, 2) << "],[" << A(2, 0)
              << "," << A(2, 1) << "," << A(2, 2) << "]]";
  }
};
template <typename T, dnup_t dnup1, dnup_t dnup2>
struct nan<mat3<T, dnup1, dnup2> > {
  constexpr mat3<T, dnup1, dnup2> operator()() const {
    return mat3<T, dnup1, dnup2>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup1, Weyl::dnup_t dnup2>
struct zero<Weyl::mat3<T, dnup1, dnup2> > {
  constexpr Weyl::mat3<T, dnup1, dnup2> operator()() const {
    return Weyl::mat3<T, dnup1, dnup2>(zero<vect<T, 6> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct norm1<mat3<T, dnup1, dnup2> > {
  typedef typename norm1<vect<T, 6> >::result_type result_type;
  constexpr result_type operator()(const mat3<T, dnup1, dnup2> &x) const {
    return norm1<vect<T, 6> >()(x.elts);
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup3>
constexpr mat3<T, dnup1, dnup2> mul(const mat3<T, dnup1, dnup3> &A,
                                    const mat3<T, !dnup3, dnup2> &B) {
  // C[a,b] = A[a,c] B[c,b]
  return mat3<T, dnup1, dnup2>([&](int a, int b) {
    return sum1([&](int x) { return A(a, x) * B(x, b); });
  });
}

template <typename F,
          typename R = remove_cv_t<remove_reference_t<result_of_t<F(int)> > > >
constexpr R sum1(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 3; ++x)
    s += f(x);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int)> > > >
constexpr R sum2(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 3; ++x)
    for (int y = 0; y < 3; ++y)
      s += f(x, y);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int, int)> > > >
constexpr R sum3(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 3; ++x)
    for (int y = 0; y < 3; ++y)
      for (int z = 0; z < 3; ++z)
        s += f(x, y, z);
  return s;
}

////////////////////////////////////////////////////////////////////////////////

// 4-vector
template <typename T, dnup_t dnup> class vec4 {
  vect<T, 4> elts;

  static constexpr int ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 4);
#endif
    return n;
  }

public:
  explicit constexpr vec4() : elts{nan<vect<T, 4> >()()} {}

  constexpr vec4(const vect<T, 4> &elts) : elts(elts) {}
  constexpr vec4(vect<T, 4> &&elts) : elts(move(elts)) {}

  explicit constexpr vec4(T vt, T vx, T vy, T vz)
      : elts(make_tuple(move(vt), move(vx), move(vy), move(vz))) {}

  constexpr vec4(initializer_list<T> v) : elts(v) {}
  // constexpr vec4(const vector<T> &v) : elts(v) {}
  // constexpr vec4(vector<T> &&v) : elts(move(v)) {}

  vec4(const GF3D2<add_const_t<T> > &gf_vt_,
       const GF3D2<add_const_t<T> > &gf_vx_,
       const GF3D2<add_const_t<T> > &gf_vy_,
       const GF3D2<add_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec4{gf_vt_(I), gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  vec4(const GF3D2<remove_const_t<T> > &gf_vt_,
       const GF3D2<remove_const_t<T> > &gf_vx_,
       const GF3D2<remove_const_t<T> > &gf_vy_,
       const GF3D2<remove_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec4{gf_vt_(I), gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  template <typename F, typename = result_of_t<F(int)> >
  constexpr vec4(const F &f) : elts{f(0), f(1), f(2), f(3)} {}

  void store(const GF3D2<T> &gf_vt_, const GF3D2<T> &gf_vx_,
             const GF3D2<T> &gf_vy_, const GF3D2<T> &gf_vz_,
             const vect<int, 3> &I) const {
    const auto &v = *this;
#ifdef CCTK_DEBUG
    if (!((CCTK_isfinite(v(0))) && (CCTK_isfinite(v(1))) &&
          (CCTK_isfinite(v(2))) && (CCTK_isfinite(v(3))))) {
      ostringstream buf;
      buf << "v=" << v;
      CCTK_VERROR("nan found: %s", buf.str().c_str());
    }
    assert(CCTK_isfinite(v(0)));
    assert(CCTK_isfinite(v(1)));
    assert(CCTK_isfinite(v(2)));
    assert(CCTK_isfinite(v(3)));
#endif
    gf_vt_(I) = v(0);
    gf_vx_(I) = v(1);
    gf_vy_(I) = v(2);
    gf_vz_(I) = v(3);
  }

  const T &operator()(int i) const { return elts[ind(i)]; }
  T &operator()(int i) { return elts[ind(i)]; }

  template <typename U = T>
  vec4<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >, dnup>
  operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I), elts[3](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  vec4<remove_cv_t<
           remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 4>)> > >,
       dnup>
  operator()(const GF3D5layout &layout, const vect<int, 4> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[3](layout, I)};
  }

  friend constexpr vec4<T, dnup> operator+(const vec4<T, dnup> &x) {
    return {+x.elts};
  }
  friend constexpr vec4<T, dnup> operator-(const vec4<T, dnup> &x) {
    return {-x.elts};
  }
  friend constexpr vec4<T, dnup> operator+(const vec4<T, dnup> &x,
                                           const vec4<T, dnup> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr vec4<T, dnup> operator-(const vec4<T, dnup> &x,
                                           const vec4<T, dnup> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr vec4<T, dnup> operator*(const T &a, const vec4<T, dnup> &x) {
    return {a * x.elts};
  }
  friend constexpr vec4<T, dnup> operator*(const vec4<T, dnup> &x, const T &a) {
    return {x.elts * a};
  }
  friend constexpr vec4<T, dnup> operator/(const vec4<T, dnup> &x, const T &a) {
    return {x.elts / a};
  }

  vec4<T, dnup> &operator*=(const T &a) { return *this = *this * a; }
  vec4<T, dnup> &operator/=(const T &a) { return *this = *this / a; }

  friend constexpr bool operator==(const vec4<T, dnup> &x,
                                   const vec4<T, dnup> &y) {
    return equal_to<vect<T, 4> >()(x.elts, y.elts);
  }
  friend constexpr bool operator!=(const vec4<T, dnup> &x,
                                   const vec4<T, dnup> &y) {
    return !(x == y);
  }

  constexpr vec4<bool, dnup> isnan1() const {
    using std::isnan1;
    return vec4<bool, dnup>(isnan1(elts));
  }

  constexpr bool any() const { return elts.any(); }
  constexpr T maxabs() const { return elts.maxabs(); }

  friend struct norm1<vec4>;

  friend ostream &operator<<(ostream &os, const vec4<T, dnup> &v) {
    return os << "[" << v(0) << "," << v(1) << "," << v(2) << "," << v(3)
              << "]";
  }
};

template <typename T, dnup_t dnup> struct nan<vec4<T, dnup> > {
  constexpr vec4<T, dnup> operator()() const { return vec4<T, dnup>(); }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup> struct zero<Weyl::vec4<T, dnup> > {
  constexpr Weyl::vec4<T, dnup> operator()() const {
    return Weyl::vec4<T, dnup>(zero<vect<T, 4> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup> struct norm1<vec4<T, dnup> > {
  typedef typename norm1<vect<T, 4> >::result_type result_type;
  constexpr result_type operator()(const vec4<T, dnup> &x) const {
    return norm1<vect<T, 4> >()(x.elts);
  }
};

////////////////////////////////////////////////////////////////////////////////

// Symmetric 4-matrix
template <typename T, dnup_t dnup1, dnup_t dnup2> class mat4 {
  static_assert(dnup1 == dnup2, "");

  vect<T, 10> elts;

  static constexpr int symind(const int i, const int j) {
#ifdef CCTK_DEBUG
    assert(i >= 0 && i <= j && j < 4);
#endif
    const int n = i * (7 - i) / 2 + j;
    // i j n
    // 0 0 0
    // 0 1 1
    // 0 2 2
    // 0 3 3
    // 1 1 4
    // 1 2 5
    // 1 3 6
    // 2 2 7
    // 2 3 8
    // 3 3 9
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 10);
#endif
    return n;
  }
  static constexpr int ind(const int i, const int j) {
    return symind(min(i, j), max(i, j));
  }

  static_assert(symind(0, 0) == 0, "");
  static_assert(symind(0, 1) == 1, "");
  static_assert(symind(0, 2) == 2, "");
  static_assert(symind(0, 3) == 3, "");
  static_assert(symind(1, 1) == 4, "");
  static_assert(symind(1, 2) == 5, "");
  static_assert(symind(1, 3) == 6, "");
  static_assert(symind(2, 2) == 7, "");
  static_assert(symind(2, 3) == 8, "");
  static_assert(symind(3, 3) == 9, "");

  // nvcc doesn't handle these constexpr expressions
#ifndef __CUDACC__
  static_assert(ind(1, 0) == ind(0, 1), "");
  static_assert(ind(2, 0) == ind(0, 2), "");
  static_assert(ind(3, 0) == ind(0, 3), "");
  static_assert(ind(2, 1) == ind(1, 2), "");
  static_assert(ind(3, 1) == ind(1, 3), "");
  static_assert(ind(3, 2) == ind(2, 3), "");
#endif

public:
  explicit constexpr mat4() : elts{nan<vect<T, 10> >()()} {}

  constexpr mat4(const vect<T, 10> &elts) : elts(elts) {}
  constexpr mat4(vect<T, 10> &&elts) : elts(move(elts)) {}

  explicit constexpr mat4(T Att, T Atx, T Aty, T Atz, T Axx, T Axy, T Axz,
                          T Ayy, T Ayz, T Azz)
      : elts(make_tuple(move(Att), move(Atx), move(Aty), move(Atz), move(Axx),
                        move(Axy), move(Axz), move(Ayy), move(Ayz),
                        move(Azz))) {}

  constexpr mat4(initializer_list<T> A) : elts(A) {}
  // constexpr mat4(const vector<T> &A) : elts(A) {}
  // constexpr mat4(vector<T> &&A) : elts(move(A)) {}

  mat4(const GF3D2<add_const_t<T> > &gf_Att_,
       const GF3D2<add_const_t<T> > &gf_Atx_,
       const GF3D2<add_const_t<T> > &gf_Aty_,
       const GF3D2<add_const_t<T> > &gf_Atz_,
       const GF3D2<add_const_t<T> > &gf_Axx_,
       const GF3D2<add_const_t<T> > &gf_Axy_,
       const GF3D2<add_const_t<T> > &gf_Axz_,
       const GF3D2<add_const_t<T> > &gf_Ayy_,
       const GF3D2<add_const_t<T> > &gf_Ayz_,
       const GF3D2<add_const_t<T> > &gf_Azz_, const vect<int, 3> &I)
      : mat4{gf_Att_(I), gf_Atx_(I), gf_Aty_(I), gf_Atz_(I), gf_Axx_(I),
             gf_Axy_(I), gf_Axz_(I), gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  mat4(const GF3D2<remove_const_t<T> > &gf_Att_,
       const GF3D2<remove_const_t<T> > &gf_Atx_,
       const GF3D2<remove_const_t<T> > &gf_Aty_,
       const GF3D2<remove_const_t<T> > &gf_Atz_,
       const GF3D2<remove_const_t<T> > &gf_Axx_,
       const GF3D2<remove_const_t<T> > &gf_Axy_,
       const GF3D2<remove_const_t<T> > &gf_Axz_,
       const GF3D2<remove_const_t<T> > &gf_Ayy_,
       const GF3D2<remove_const_t<T> > &gf_Ayz_,
       const GF3D2<remove_const_t<T> > &gf_Azz_, const vect<int, 3> &I)
      : mat4{gf_Att_(I), gf_Atx_(I), gf_Aty_(I), gf_Atz_(I), gf_Axx_(I),
             gf_Axy_(I), gf_Axz_(I), gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  template <typename F, typename = result_of_t<F(int, int)> >
  constexpr mat4(const F &f)
      : elts{f(0, 0), f(0, 1), f(0, 2), f(0, 3), f(1, 1),
             f(1, 2), f(1, 3), f(2, 2), f(2, 3), f(3, 3)} {
#ifdef CCTK_DEBUG
    // Check symmetry
    const T f10 = f(1, 0);
    const T f20 = f(2, 0);
    const T f30 = f(3, 0);
    const T f21 = f(2, 1);
    const T f31 = f(3, 1);
    const T f32 = f(3, 2);
    const auto scale = norm1<T>()(elts[0]) + norm1<T>()(elts[1]) +
                       norm1<T>()(elts[2]) + norm1<T>()(elts[3]) +
                       norm1<T>()(elts[4]) + norm1<T>()(elts[5]) +
                       norm1<T>()(elts[6]) + norm1<T>()(elts[7]) +
                       norm1<T>()(elts[8]) + norm1<T>()(elts[9]);
    const auto is_approx{[scale](const T &fgood, const T &fother) {
      return norm1<T>()(fother - fgood) <= 1.0e-12 * (1 + scale);
    }};
    if (!(is_approx(f10, (*this)(0, 1)) && is_approx(f20, (*this)(0, 2)) &&
          is_approx(f30, (*this)(0, 3)) && is_approx(f21, (*this)(1, 2)) &&
          is_approx(f31, (*this)(1, 3)) && is_approx(f32, (*this)(2, 3)))) {
      ostringstream buf;
      buf << "f(0,1)=" << (*this)(0, 1) << "\n"
          << "f(1,0)=" << f10 << "\n"
          << "f(0,2)=" << (*this)(0, 2) << "\n"
          << "f(2,0)=" << f20 << "\n"
          << "f(0,3)=" << (*this)(0, 3) << "\n"
          << "f(3,0)=" << f30 << "\n"
          << "f(1,2)=" << (*this)(1, 2) << "\n"
          << "f(2,1)=" << f21 << "\n"
          << "f(1,3)=" << (*this)(1, 3) << "\n"
          << "f(3,1)=" << f31 << "\n"
          << "f(2,3)=" << (*this)(2, 3) << "\n"
          << "f(3,2)=" << f32 << "\n";
      CCTK_VERROR("symmetric matrix is not symmetric:\n%s", buf.str().c_str());
    }
    assert(is_approx((*this)(0, 1), f10));
    assert(is_approx((*this)(0, 2), f20));
    assert(is_approx((*this)(0, 3), f30));
    assert(is_approx((*this)(1, 2), f21));
    assert(is_approx((*this)(1, 3), f31));
    assert(is_approx((*this)(2, 3), f32));
#endif
  }

  void store(const GF3D2<T> &gf_Att_, const GF3D2<T> &gf_Atx_,
             const GF3D2<T> &gf_Aty_, const GF3D2<T> &gf_Atz_,
             const GF3D2<T> &gf_Axx_, const GF3D2<T> &gf_Axy_,
             const GF3D2<T> &gf_Axz_, const GF3D2<T> &gf_Ayy_,
             const GF3D2<T> &gf_Ayz_, const GF3D2<T> &gf_Azz_,
             const vect<int, 3> &I) const {
    const auto &A = *this;
#ifdef CCTK_DEBUG
    assert(CCTK_isfinite(A(0, 0)));
    assert(CCTK_isfinite(A(0, 1)));
    assert(CCTK_isfinite(A(0, 2)));
    assert(CCTK_isfinite(A(0, 3)));
    assert(CCTK_isfinite(A(1, 1)));
    assert(CCTK_isfinite(A(1, 2)));
    assert(CCTK_isfinite(A(1, 3)));
    assert(CCTK_isfinite(A(2, 2)));
    assert(CCTK_isfinite(A(2, 3)));
    assert(CCTK_isfinite(A(3, 3)));
#endif
    gf_Att_(I) = A(0, 0);
    gf_Atx_(I) = A(0, 1);
    gf_Aty_(I) = A(0, 2);
    gf_Atz_(I) = A(0, 3);
    gf_Axx_(I) = A(1, 1);
    gf_Axy_(I) = A(1, 2);
    gf_Axz_(I) = A(1, 3);
    gf_Ayy_(I) = A(2, 2);
    gf_Ayz_(I) = A(2, 3);
    gf_Azz_(I) = A(3, 3);
  }

  const T &operator()(int i, int j) const { return elts[ind(i, j)]; }
  //  T &operator()(int i, int j) { return
  // elts[symind(i, j)]; }
  T &operator()(int i, int j) { return elts[ind(i, j)]; }

  template <typename U = T>
  mat4<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 10>)> > >, dnup1,
       dnup2>
  operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I), elts[4](I), elts[4](I),
            elts[5](I), elts[6](I), elts[7](I), elts[8](I), elts[9](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  mat4<remove_cv_t<
           remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
       dnup1, dnup2>
  operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[4](layout, I), elts[4](layout, I), elts[5](layout, I),
            elts[6](layout, I), elts[7](layout, I), elts[8](layout, I),
            elts[9](layout, I)};
  }

  friend constexpr mat4<T, dnup1, dnup2>
  operator+(const mat4<T, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend constexpr mat4<T, dnup1, dnup2>
  operator-(const mat4<T, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend constexpr mat4<T, dnup1, dnup2>
  operator+(const mat4<T, dnup1, dnup2> &x, const mat4<T, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr mat4<T, dnup1, dnup2>
  operator-(const mat4<T, dnup1, dnup2> &x, const mat4<T, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr mat4<T, dnup1, dnup2>
  operator*(const T &a, const mat4<T, dnup1, dnup2> &x) {
    return {a * x.elts};
  }

  friend constexpr bool operator==(const mat4<T, dnup1, dnup2> &x,
                                   const mat4<T, dnup1, dnup2> &y) {
    return equal_to<vect<T, 10> >()(x.elts, y.elts);
  }
  friend constexpr bool operator!=(const mat4<T, dnup1, dnup2> &x,
                                   const mat4<T, dnup1, dnup2> &y) {
    return !(x == y);
  }

  constexpr T maxabs() const { return elts.maxabs(); }

  friend struct norm1<mat4>;

  constexpr T det() const {
    const auto &A = *this;
    return -(A(0, 0) * pow(A(1, 3), 2) * A(2, 2)) +
           pow(A(0, 3), 2) * (pow(A(1, 2), 2) - A(1, 1) * A(2, 2)) +
           2 * A(0, 0) * A(1, 2) * A(1, 3) * A(2, 3) +
           pow(A(0, 1), 2) * pow(A(2, 3), 2) -
           A(0, 0) * A(1, 1) * pow(A(2, 3), 2) -
           2 * A(0, 3) *
               (A(0, 2) * (A(1, 2) * A(1, 3) - A(1, 1) * A(2, 3)) +
                A(0, 1) * (-(A(1, 3) * A(2, 2)) + A(1, 2) * A(2, 3))) -
           A(0, 0) * pow(A(1, 2), 2) * A(3, 3) -
           pow(A(0, 1), 2) * A(2, 2) * A(3, 3) +
           A(0, 0) * A(1, 1) * A(2, 2) * A(3, 3) +
           pow(A(0, 2), 2) * (pow(A(1, 3), 2) - A(1, 1) * A(3, 3)) +
           2 * A(0, 1) * A(0, 2) * (-(A(1, 3) * A(2, 3)) + A(1, 2) * A(3, 3));
  }

  constexpr mat4<T, !dnup1, !dnup2> inv(const T detA) const {
    const auto &A = *this;
    const T detA1 = 1 / detA;
    return mat4<T, !dnup1, !dnup2>{
        detA1 * (-(pow(A(1, 3), 2) * A(2, 2)) +
                 2 * A(1, 2) * A(1, 3) * A(2, 3) - pow(A(1, 2), 2) * A(3, 3) +
                 A(1, 1) * (-pow(A(2, 3), 2) + A(2, 2) * A(3, 3))),
        detA1 * (A(0, 3) * (A(1, 3) * A(2, 2) - A(1, 2) * A(2, 3)) +
                 A(0, 2) * (-(A(1, 3) * A(2, 3)) + A(1, 2) * A(3, 3)) +
                 A(0, 1) * (pow(A(2, 3), 2) - A(2, 2) * A(3, 3))),
        detA1 * (A(0, 3) * (-(A(1, 2) * A(1, 3)) + A(1, 1) * A(2, 3)) +
                 A(0, 2) * (pow(A(1, 3), 2) - A(1, 1) * A(3, 3)) +
                 A(0, 1) * (-(A(1, 3) * A(2, 3)) + A(1, 2) * A(3, 3))),
        detA1 * (A(0, 3) * (pow(A(1, 2), 2) - A(1, 1) * A(2, 2)) +
                 A(0, 2) * (-(A(1, 2) * A(1, 3)) + A(1, 1) * A(2, 3)) +
                 A(0, 1) * (A(1, 3) * A(2, 2) - A(1, 2) * A(2, 3))),
        detA1 * (-(pow(A(0, 3), 2) * A(2, 2)) +
                 2 * A(0, 2) * A(0, 3) * A(2, 3) - pow(A(0, 2), 2) * A(3, 3) +
                 A(0, 0) * (-pow(A(2, 3), 2) + A(2, 2) * A(3, 3))),
        detA1 * (pow(A(0, 3), 2) * A(1, 2) -
                 A(0, 3) * (A(0, 2) * A(1, 3) + A(0, 1) * A(2, 3)) +
                 A(0, 1) * A(0, 2) * A(3, 3) +
                 A(0, 0) * (A(1, 3) * A(2, 3) - A(1, 2) * A(3, 3))),
        detA1 * (pow(A(0, 2), 2) * A(1, 3) + A(0, 1) * A(0, 3) * A(2, 2) -
                 A(0, 2) * (A(0, 3) * A(1, 2) + A(0, 1) * A(2, 3)) +
                 A(0, 0) * (-(A(1, 3) * A(2, 2)) + A(1, 2) * A(2, 3))),
        detA1 * (-(pow(A(0, 3), 2) * A(1, 1)) +
                 2 * A(0, 1) * A(0, 3) * A(1, 3) - pow(A(0, 1), 2) * A(3, 3) +
                 A(0, 0) * (-pow(A(1, 3), 2) + A(1, 1) * A(3, 3))),
        detA1 * (-(A(0, 1) * A(0, 3) * A(1, 2)) +
                 A(0, 2) * (A(0, 3) * A(1, 1) - A(0, 1) * A(1, 3)) +
                 pow(A(0, 1), 2) * A(2, 3) +
                 A(0, 0) * (A(1, 2) * A(1, 3) - A(1, 1) * A(2, 3))),
        detA1 * (-(pow(A(0, 2), 2) * A(1, 1)) +
                 2 * A(0, 1) * A(0, 2) * A(1, 2) - pow(A(0, 1), 2) * A(2, 2) +
                 A(0, 0) * (-pow(A(1, 2), 2) + A(1, 1) * A(2, 2)))};

    return mat4<T, !dnup1, !dnup2>{
        detA1 * (-(pow(A(1, 3), 2) * A(2, 2)) +
                 2 * A(1, 2) * A(1, 3) * A(2, 3) - pow(A(1, 2), 2) * A(3, 3) +
                 A(1, 1) * (-pow(A(2, 3), 2) + A(2, 2) * A(3, 3))),
        detA1 * (A(0, 3) * (A(1, 3) * A(2, 2) - A(1, 2) * A(2, 3)) +
                 A(0, 2) * (-(A(1, 3) * A(2, 3)) + A(1, 2) * A(3, 3)) +
                 A(0, 1) * (pow(A(2, 3), 2) - A(2, 2) * A(3, 3))),
        detA1 * (A(0, 3) * (-(A(1, 2) * A(1, 3)) + A(1, 1) * A(2, 3)) +
                 A(0, 2) * (pow(A(1, 3), 2) - A(1, 1) * A(3, 3)) +
                 A(0, 1) * (-(A(1, 3) * A(2, 3)) + A(1, 2) * A(3, 3))),
        detA1 * (A(0, 3) * (pow(A(1, 2), 2) - A(1, 1) * A(2, 2)) +
                 A(0, 2) * (-(A(1, 2) * A(1, 3)) + A(1, 1) * A(2, 3)) +
                 A(0, 1) * (A(1, 3) * A(2, 2) - A(1, 2) * A(2, 3))),
        detA1 * (-(pow(A(0, 3), 2) * A(2, 2)) +
                 2 * A(0, 2) * A(0, 3) * A(2, 3) - pow(A(0, 2), 2) * A(3, 3) +
                 A(0, 0) * (-pow(A(2, 3), 2) + A(2, 2) * A(3, 3))),
        detA1 * (pow(A(0, 3), 2) * A(1, 2) -
                 A(0, 3) * (A(0, 2) * A(1, 3) + A(0, 1) * A(2, 3)) +
                 A(0, 1) * A(0, 2) * A(3, 3) +
                 A(0, 0) * (A(1, 3) * A(2, 3) - A(1, 2) * A(3, 3))),
        detA1 * (pow(A(0, 2), 2) * A(1, 3) + A(0, 1) * A(0, 3) * A(2, 2) -
                 A(0, 2) * (A(0, 3) * A(1, 2) + A(0, 1) * A(2, 3)) +
                 A(0, 0) * (-(A(1, 3) * A(2, 2)) + A(1, 2) * A(2, 3))),
        detA1 * (-(pow(A(0, 3), 2) * A(1, 1)) +
                 2 * A(0, 1) * A(0, 3) * A(1, 3) - pow(A(0, 1), 2) * A(3, 3) +
                 A(0, 0) * (-pow(A(1, 3), 2) + A(1, 1) * A(3, 3))),
        detA1 * (-(A(0, 1) * A(0, 3) * A(1, 2)) +
                 A(0, 2) * (A(0, 3) * A(1, 1) - A(0, 1) * A(1, 3)) +
                 pow(A(0, 1), 2) * A(2, 3) +
                 A(0, 0) * (A(1, 2) * A(1, 3) - A(1, 1) * A(2, 3))),
        detA1 * (-(pow(A(0, 2), 2) * A(1, 1)) +
                 2 * A(0, 1) * A(0, 2) * A(1, 2) - pow(A(0, 1), 2) * A(2, 2) +
                 A(0, 0) * (-pow(A(1, 2), 2) + A(1, 1) * A(2, 2)))};
  }

  constexpr T trace(const mat4<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    return sum2([&](int x, int y) { return gu(x, y) * A(x, y); });
  }

#if 0
  constexpr  mat4 trace_free(
      const mat4<T, dnup1, dnup2> &g, const mat4<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    const T trA = A.trace(gu);
    return mat4([&](int a, int b)  {
      return A(a, b) - trA / 2 * g(a, b);
    });
  }
#endif

  friend ostream &operator<<(ostream &os, const mat4<T, dnup1, dnup2> &A) {
    return os << "[[" << A(0, 0) << "," << A(0, 1) << "," << A(0, 2) << ","
              << A(0, 3) << "],[" << A(1, 0) << "," << A(1, 1) << "," << A(1, 2)
              << "," << A(1, 3) << "],[" << A(2, 0) << "," << A(2, 1) << ","
              << A(2, 2) << "," << A(2, 3) << "],[" << A(3, 0) << "," << A(3, 1)
              << "," << A(3, 2) << "," << A(3, 3) << "]]";
  }
};
template <typename T, dnup_t dnup1, dnup_t dnup2>
struct nan<mat4<T, dnup1, dnup2> > {
  constexpr mat4<T, dnup1, dnup2> operator()() const {
    return mat4<T, dnup1, dnup2>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup1, Weyl::dnup_t dnup2>
struct zero<Weyl::mat4<T, dnup1, dnup2> > {
  constexpr Weyl::mat4<T, dnup1, dnup2> operator()() const {
    return Weyl::mat4<T, dnup1, dnup2>(zero<vect<T, 10> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct norm1<mat4<T, dnup1, dnup2> > {
  typedef typename norm1<vect<T, 10> >::result_type result_type;
  constexpr result_type operator()(const mat4<T, dnup1, dnup2> &x) const {
    return norm1<vect<T, 10> >()(x.elts);
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup4>
constexpr mat4<T, dnup1, dnup2> mul(const mat4<T, dnup1, dnup4> &A,
                                    const mat4<T, !dnup4, dnup2> &B) {
  // C[a,b] = A[a,c] B[c,b]
  return mat4<T, dnup1, dnup2>([&](int a, int b) {
    return sum1([&](int x) { return A(a, x) * B(x, b); });
  });
}

////////////////////////////////////////////////////////////////////////////////

// Antisymmetric 4-matrix
template <typename T, dnup_t dnup1, dnup_t dnup2> class amat4 {
  static_assert(dnup1 == dnup2, "");

  vect<T, 6> elts;

  static constexpr int asymind(const int i, const int j) {
#ifdef CCTK_DEBUG
    assert(i >= 0 && i < j && j < 4);
#endif
    const int n = i * (5 - i) / 2 + j - 1;
    // i j n
    // 0 1 0
    // 0 2 1
    // 0 3 2
    // 1 2 3
    // 1 3 4
    // 2 3 5
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 6);
#endif
    return n;
  }
  static constexpr int ind(const int i, const int j) {
    return asymind(min(i, j), max(i, j));
  }
  static constexpr int sign(const int i, const int j) {
    if (i == j)
      return 0;
    return i < j ? 1 : -1;
  }

  static_assert(asymind(0, 1) == 0, "");
  static_assert(asymind(0, 2) == 1, "");
  static_assert(asymind(0, 3) == 2, "");
  static_assert(asymind(1, 2) == 3, "");
  static_assert(asymind(1, 3) == 4, "");
  static_assert(asymind(2, 3) == 5, "");

  // nvcc doesn't handle these constexpr expressions
#ifndef __CUDACC__
  static_assert(ind(1, 0) == ind(0, 1), "");
  static_assert(ind(2, 0) == ind(0, 2), "");
  static_assert(ind(3, 0) == ind(0, 3), "");
  static_assert(ind(2, 1) == ind(1, 2), "");
  static_assert(ind(3, 1) == ind(1, 3), "");
  static_assert(ind(3, 2) == ind(2, 3), "");
#endif

  static_assert(sign(1, 0) == -sign(0, 1), "");
  static_assert(sign(2, 0) == -sign(0, 2), "");
  static_assert(sign(3, 0) == -sign(0, 3), "");
  static_assert(sign(2, 1) == -sign(1, 2), "");
  static_assert(sign(3, 1) == -sign(1, 3), "");
  static_assert(sign(3, 2) == -sign(2, 3), "");
  static_assert(sign(0, 0) == 0, "");
  static_assert(sign(1, 1) == 0, "");
  static_assert(sign(2, 2) == 0, "");
  static_assert(sign(3, 3) == 0, "");

public:
  explicit constexpr amat4() : elts{nan<vect<T, 6> >()()} {}

  constexpr amat4(const vect<T, 6> &elts) : elts(elts) {}
  constexpr amat4(vect<T, 6> &&elts) : elts(move(elts)) {}

  explicit constexpr amat4(T Atx, T Aty, T Atz, T Axy, T Axz, T Ayz)
      : elts(make_tuple(move(Atx), move(Aty), move(Atz), move(Axy), move(Axz),
                        move(Ayz))) {}

  constexpr amat4(initializer_list<T> A) : elts(A) {}
  // constexpr amat4(const vector<T> &A) : elts(A) {}
  // constexpr amat4(vector<T> &&A) : elts(move(A)) {}

  amat4(const GF3D2<add_const_t<T> > &gf_Atx_,
        const GF3D2<add_const_t<T> > &gf_Aty_,
        const GF3D2<add_const_t<T> > &gf_Atz_,
        const GF3D2<add_const_t<T> > &gf_Axy_,
        const GF3D2<add_const_t<T> > &gf_Axz_,
        const GF3D2<add_const_t<T> > &gf_Ayz_, const vect<int, 3> &I)
      : amat4{gf_Atx_(I), gf_Aty_(I), gf_Atz_(I),
              gf_Axy_(I), gf_Axz_(I), gf_Ayz_(I)} {}

  amat4(const GF3D2<remove_const_t<T> > &gf_Atx_,
        const GF3D2<remove_const_t<T> > &gf_Aty_,
        const GF3D2<remove_const_t<T> > &gf_Atz_,
        const GF3D2<remove_const_t<T> > &gf_Axy_,
        const GF3D2<remove_const_t<T> > &gf_Axz_,
        const GF3D2<remove_const_t<T> > &gf_Ayz_, const vect<int, 3> &I)
      : amat4{gf_Atx_(I), gf_Aty_(I), gf_Atz_(I),
              gf_Axy_(I), gf_Axz_(I), gf_Ayz_(I)} {}

  template <typename F, typename = result_of_t<F(int, int)> >
  constexpr amat4(const F &f)
      : elts{f(0, 1), f(0, 2), f(0, 3), f(1, 2), f(1, 3), f(2, 3)} {
#ifdef CCTK_DEBUG
    // Check antisymmetry
    const T f00 = f(0, 0);
    const T f10 = f(1, 0);
    const T f20 = f(2, 0);
    const T f30 = f(3, 0);
    const T f11 = f(1, 1);
    const T f21 = f(2, 1);
    const T f31 = f(3, 1);
    const T f22 = f(2, 2);
    const T f32 = f(3, 2);
    const T f33 = f(3, 3);
    const auto scale = norm1<T>()(elts[0]) + norm1<T>()(elts[1]) +
                       norm1<T>()(elts[2]) + norm1<T>()(elts[3]) +
                       norm1<T>()(elts[4]) + norm1<T>()(elts[5]);
    const auto is_approx{[scale](const T &fgood, const T &fother) {
      return norm1<T>()(fother - fgood) <= 1.0e-12 * (1 + scale);
    }};
    if (!(is_approx(f00, zero<T>()()) && is_approx(f10, -(*this)(0, 1)) &&
          is_approx(f20, -(*this)(0, 2)) && is_approx(f30, -(*this)(0, 3)) &&
          is_approx(f11, zero<T>()()) && is_approx(f21, -(*this)(1, 2)) &&
          is_approx(f31, -(*this)(1, 3)) && is_approx(f22, zero<T>()()) &&
          is_approx(f32, -(*this)(2, 3)) && is_approx(f33, zero<T>()()))) {
      ostringstream buf;
      buf << "f(0,0)=" << f00 << "\n"
          << "f(0,1)=" << (*this)(0, 1) << "\n"
          << "f(1,0)=" << f10 << "\n"
          << "f(0,2)=" << (*this)(0, 2) << "\n"
          << "f(2,0)=" << f20 << "\n"
          << "f(0,3)=" << (*this)(0, 3) << "\n"
          << "f(3,0)=" << f30 << "\n"
          << "f(1,1)=" << f11 << "\n"
          << "f(1,2)=" << (*this)(1, 2) << "\n"
          << "f(2,1)=" << f21 << "\n"
          << "f(1,3)=" << (*this)(1, 3) << "\n"
          << "f(3,1)=" << f31 << "\n"
          << "f(2,2)=" << f22 << "\n"
          << "f(2,3)=" << (*this)(2, 3) << "\n"
          << "f(3,2)=" << f32 << "\n"
          << "f(3,3)=" << f33 << "\n";
      CCTK_VERROR("antisymmetric matrix is not antisymmetric:\n%s",
                  buf.str().c_str());
    }
    assert(is_approx(f00, zero<T>()()));
    assert(is_approx(f10, -(*this)(0, 1)));
    assert(is_approx(f20, -(*this)(0, 2)));
    assert(is_approx(f30, -(*this)(0, 3)));
    assert(is_approx(f11, zero<T>()()));
    assert(is_approx(f21, -(*this)(1, 2)));
    assert(is_approx(f31, -(*this)(1, 3)));
    assert(is_approx(f22, zero<T>()()));
    assert(is_approx(f32, -(*this)(2, 3)));
    assert(is_approx(f33, zero<T>()()));
#endif
  }

  void store(const GF3D2<T> &gf_Atx_, const GF3D2<T> &gf_Aty_,
             const GF3D2<T> &gf_Atz_, const GF3D2<T> &gf_Axy_,
             const GF3D2<T> &gf_Axz_, const GF3D2<T> &gf_Ayz_,
             const vect<int, 3> &I) const {
    const auto &A = *this;
#ifdef CCTK_DEBUG
    assert(CCTK_isfinite(A(0, 1)));
    assert(CCTK_isfinite(A(0, 2)));
    assert(CCTK_isfinite(A(0, 3)));
    assert(CCTK_isfinite(A(1, 2)));
    assert(CCTK_isfinite(A(1, 3)));
    assert(CCTK_isfinite(A(2, 3)));
#endif
    gf_Atx_(I) = A(0, 1);
    gf_Aty_(I) = A(0, 2);
    gf_Atz_(I) = A(0, 3);
    gf_Axy_(I) = A(1, 2);
    gf_Axz_(I) = A(1, 3);
    gf_Ayz_(I) = A(2, 3);
  }

  const T operator()(int i, int j) const {
    const int s = sign(i, j);
    if (s == 0)
      return zero<T>()();
    return s * elts[ind(i, j)];
  }
  //  T &operator()(int i, int j) { return
  // elts[symind(i, j)]; }
  T &operator()(int i, int j) {
    const int s = sign(i, j);
    assert(s == 1);
    return elts[ind(i, j)];
  }

  template <typename U = T>
  amat4<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 6>)> > >, dnup1,
        dnup2>
  operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I),
            elts[3](I), elts[4](I), elts[5](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  amat4<remove_cv_t<
            remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
        dnup1, dnup2>
  operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[4](layout, I), elts[4](layout, I), elts[5](layout, I)};
  }

  friend constexpr amat4<T, dnup1, dnup2>
  operator+(const amat4<T, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend constexpr amat4<T, dnup1, dnup2>
  operator-(const amat4<T, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend constexpr amat4<T, dnup1, dnup2>
  operator+(const amat4<T, dnup1, dnup2> &x, const amat4<T, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr amat4<T, dnup1, dnup2>
  operator-(const amat4<T, dnup1, dnup2> &x, const amat4<T, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr amat4<T, dnup1, dnup2>
  operator*(const T &a, const amat4<T, dnup1, dnup2> &x) {
    return {a * x.elts};
  }

  friend constexpr bool operator==(const amat4<T, dnup1, dnup2> &x,
                                   const amat4<T, dnup1, dnup2> &y) {
    return equal_to<vect<T, 6> >()(x.elts, y.elts);
  }
  friend constexpr bool operator!=(const amat4<T, dnup1, dnup2> &x,
                                   const amat4<T, dnup1, dnup2> &y) {
    return !(x == y);
  }

  constexpr T maxabs() const { return elts.maxabs(); }

  friend struct norm1<amat4>;

  constexpr T det() const {
    const auto &A = *this;
    return pow(A(0, 3), 2) * pow(A(1, 2), 2) -
           2 * A(0, 2) * A(0, 3) * A(1, 2) * A(1, 3) +
           pow(A(0, 2), 2) * pow(A(1, 3), 2) +
           2 * A(0, 1) * A(0, 3) * A(1, 2) * A(2, 3) -
           2 * A(0, 1) * A(0, 2) * A(1, 3) * A(2, 3) +
           pow(A(0, 1), 2) * pow(A(2, 3), 2);
  }

  constexpr amat4<T, !dnup1, !dnup2> inv(const T detA) const {
    const auto &A = *this;
    const T detA1 = 1 / detA;
    return amat4<T, !dnup1, !dnup2>{
        detA1 * (-(A(2, 3) * (A(0, 3) * A(1, 2) - A(0, 2) * A(1, 3) +
                              A(0, 1) * A(2, 3)))),
        detA1 * (A(1, 3) *
                 (A(0, 3) * A(1, 2) - A(0, 2) * A(1, 3) + A(0, 1) * A(2, 3))),
        detA1 * (-(A(1, 2) * (A(0, 3) * A(1, 2) - A(0, 2) * A(1, 3) +
                              A(0, 1) * A(2, 3)))),
        detA1 * (-(A(0, 3) * (A(0, 3) * A(1, 2) - A(0, 2) * A(1, 3) +
                              A(0, 1) * A(2, 3)))),
        detA1 * (A(0, 2) *
                 (A(0, 3) * A(1, 2) - A(0, 2) * A(1, 3) + A(0, 1) * A(2, 3))),
        detA1 * (-(A(0, 1) * (A(0, 3) * A(1, 2) - A(0, 2) * A(1, 3) +
                              A(0, 1) * A(2, 3))))};
  }

  constexpr T trace(const amat4<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    return sum2([&](int x, int y) { return gu(x, y) * A(x, y); });
  }

#if 0
  constexpr  amat4 trace_free(
      const amat4<T, dnup1, dnup2> &g, const amat4<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    const T trA = A.trace(gu);
    return amat4([&](int a, int b)  {
      return A(a, b) - trA / 2 * g(a, b);
    });
  }
#endif

  friend ostream &operator<<(ostream &os, const amat4<T, dnup1, dnup2> &A) {
    return os << "[[" << A(0, 0) << "," << A(0, 1) << "," << A(0, 2) << ","
              << A(0, 3) << "],[" << A(1, 0) << "," << A(1, 1) << "," << A(1, 2)
              << "," << A(1, 3) << "],[" << A(2, 0) << "," << A(2, 1) << ","
              << A(2, 2) << "," << A(2, 3) << "],[" << A(3, 0) << "," << A(3, 1)
              << "," << A(3, 2) << "," << A(3, 3) << "]]";
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct nan<amat4<T, dnup1, dnup2> > {
  constexpr amat4<T, dnup1, dnup2> operator()() const {
    return amat4<T, dnup1, dnup2>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup1, Weyl::dnup_t dnup2>
struct zero<Weyl::amat4<T, dnup1, dnup2> > {
  constexpr Weyl::amat4<T, dnup1, dnup2> operator()() const {
    return Weyl::amat4<T, dnup1, dnup2>(zero<vect<T, 6> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct norm1<amat4<T, dnup1, dnup2> > {
  typedef typename norm1<vect<T, 6> >::result_type result_type;
  constexpr result_type operator()(const amat4<T, dnup1, dnup2> &x) const {
    return norm1<vect<T, 6> >()(x.elts);
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup4>
constexpr amat4<T, dnup1, dnup2> mul(const amat4<T, dnup1, dnup4> &A,
                                     const amat4<T, !dnup4, dnup2> &B) {
  // C[a,b] = A[a,c] B[c,b]
  return amat4<T, dnup1, dnup2>([&](int a, int b) {
    return sum1([&](int x) { return A(a, x) * B(x, b); });
  });
}

} // namespace Weyl

#endif // #ifndef TENSOR_HXX
