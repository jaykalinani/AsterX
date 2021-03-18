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

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(const T x) {
  return x * x;
}
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow3(const T x) {
  const T x2 = x * x;
  return x2 * x;
}
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow4(const T x) {
  const T x2 = x * x;
  return x2 * x2;
}
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow5(const T x) {
  const T x2 = x * x;
  const T x4 = x2 * x2;
  return x4 * x;
}
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow6(const T x) {
  const T x2 = x * x;
  const T x4 = x2 * x2;
  return x4 * x2;
}

namespace detail {
template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pown(const T x,
                                                                    int n) {
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

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
pown(const T x, const int n) {
  return n >= 0 ? detail::pown(x, n) : 1 / detail::pown(x, -n);
}

constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
factorial(int n) {
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
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
sum41(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 4; ++x)
    s += f(x);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int)> > > >
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
sum42(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 4; ++y)
      s += f(x, y);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int, int)> > > >
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
sum43(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 4; ++y)
      for (int z = 0; z < 4; ++z)
        s += f(x, y, z);
  return s;
}

template <typename F, typename R = remove_cv_t<remove_reference_t<
                          result_of_t<F(int, int, int, int)> > > >
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
sum44(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 4; ++x)
    for (int y = 0; y < 4; ++y)
      for (int z = 0; z < 4; ++z)
        for (int w = 0; w < 4; ++w)
          s += f(x, y, z, w);
  return s;
}

////////////////////////////////////////////////////////////////////////////////

template <typename F, typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
fold(const F &f, const T &x) {
  return x;
}
template <typename F, typename T, typename... Ts>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
fold(const F &f, const T &x0, const T &x1, const Ts &...xs) {
  return fold(f, fold(f, x0, x1), xs...);
}

template <typename T>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T add() {
  return T(0);
}
template <typename T, typename... Ts>
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
add(const T &x, const Ts &...xs) {
  return x + add(xs...);
}

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct nan {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  operator()() const {
    return NAN;
  }
};
template <typename T, int D> struct nan<vect<T, D> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vect<T, D>
  operator()() const {
    return vect<T, D>::pure(nan<T>()());
  }
};

template <typename T> struct norm1 {
  typedef T result_type;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST result_type
  operator()(const T &x) const {
    return abs(x);
  }
};
template <typename T, int D> struct norm1<vect<T, D> > {
  typedef typename norm1<T>::result_type result_type;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST result_type
  operator()(const vect<T, D> &xs) const {
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

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 3);
#endif
    return n;
  }

public:
  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3()
      : elts{nan<vect<T, 3> >()()} {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec3(vect<T, 3> elts)
      : elts(move(elts)) {}

  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec3(T vx, T vy, T vz)
      : elts(make_tuple(move(vx), move(vy), move(vz))) {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec3(initializer_list<T> v)
      : elts(v) {}
  // constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3(const
  // vector<T> &v) : elts(v) {}
  // constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  // vec3(vector<T> &&v) : elts(move(v)) {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec3(const GF3D2<add_const_t<T> > &gf_vx_,
       const GF3D2<add_const_t<T> > &gf_vy_,
       const GF3D2<add_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec3{gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec3(const GF3D2<remove_const_t<T> > &gf_vx_,
       const GF3D2<remove_const_t<T> > &gf_vy_,
       const GF3D2<remove_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec3{gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  template <typename F, typename = result_of_t<F(int)> >
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE /*CCTK_DEVICE CCTK_HOST*/
  vec3(const F &f)
      : elts{f(0), f(1), f(2)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void
  store(const GF3D2<T> &gf_vx_, const GF3D2<T> &gf_vy_, const GF3D2<T> &gf_vz_,
        const vect<int, 3> &I) const {
    const auto &v = *this;
    gf_vx_(I) = v(0);
    gf_vy_(I) = v(1);
    gf_vz_(I) = v(2);
  }

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST const T &
  operator()(int i) const {
    return elts[ind(i)];
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T &
  operator()(int i) {
    return elts[ind(i)];
  }

  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3<
      remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >, dnup>
  operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      vec3<remove_cv_t<
               remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
           dnup>
      operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I)};
  }

  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3<T, dnup>
      operator+(const vec3<T, dnup> &x) {
    return {+x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3<T, dnup>
      operator-(const vec3<T, dnup> &x) {
    return {-x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3<T, dnup>
      operator+(const vec3<T, dnup> &x, const vec3<T, dnup> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3<T, dnup>
      operator-(const vec3<T, dnup> &x, const vec3<T, dnup> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3<T, dnup>
      operator*(const T &a, const vec3<T, dnup> &x) {
    return {a * x.elts};
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator==(const vec3<T, dnup> &x, const vec3<T, dnup> &y) {
    return equal_to<vect<T, 3> >()(x.elts, y.elts);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator!=(const vec3<T, dnup> &x, const vec3<T, dnup> &y) {
    return !(x == y);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  maxabs() const {
    return elts.maxabs();
  }

  friend struct norm1<vec3>;

  friend ostream &operator<<(ostream &os, const vec3<T, dnup> &v) {
    return os << "[" << v(0) << "," << v(1) << "," << v(2) << "]";
  }
};

template <typename T, dnup_t dnup> struct nan<vec3<T, dnup> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec3<T, dnup>
  operator()() const {
    return vec3<T, dnup>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup> struct zero<Weyl::vec3<T, dnup> > {
  constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST Weyl::vec3<T, dnup>
      operator()() const {
    return Weyl::vec3<T, dnup>(zero<vect<T, 3> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup> struct norm1<vec3<T, dnup> > {
  typedef typename norm1<vect<T, 3> >::result_type result_type;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST result_type
  operator()(const vec3<T, dnup> &x) const {
    return norm1<vect<T, 3> >()(x.elts);
  }
};

////////////////////////////////////////////////////////////////////////////////

// Symmetric 3-matrix
template <typename T, dnup_t dnup1, dnup_t dnup2> class mat3 {
  static_assert(dnup1 == dnup2, "");

  vect<T, 6> elts;

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  symind(const int i, const int j) {
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
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  ind(const int i, const int j) {
    return symind(min(i, j), max(i, j));
  }

  static_assert(symind(0, 0) == 0, "");
  static_assert(symind(0, 1) == 1, "");
  static_assert(symind(0, 2) == 2, "");
  static_assert(symind(1, 1) == 3, "");
  static_assert(symind(1, 2) == 4, "");
  static_assert(symind(2, 2) == 5, "");

  // nvcc doesn't handle these constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST expressions
#ifndef __CUDACC__
  static_assert(ind(1, 0) == ind(0, 1), "");
  static_assert(ind(2, 0) == ind(0, 2), "");
  static_assert(ind(2, 1) == ind(1, 2), "");
#endif

public:
  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3()
      : elts{nan<vect<T, 6> >()()} {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat3(const vect<T, 6> &elts)
      : elts(elts) {}
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat3(vect<T, 6> &&elts)
      : elts(move(elts)) {}

  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat3(T Axx, T Axy, T Axz, T Ayy, T Ayz, T Azz)
      : elts(make_tuple(move(Axx), move(Axy), move(Axz), move(Ayy), move(Ayz),
                        move(Azz))) {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat3(initializer_list<T> A)
      : elts(A) {}
  // constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3(const
  // vector<T> &A) : elts(A) {} constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST mat3(vector<T> &&A) : elts(move(A)) {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat3(const GF3D2<add_const_t<T> > &gf_Axx_,
       const GF3D2<add_const_t<T> > &gf_Axy_,
       const GF3D2<add_const_t<T> > &gf_Axz_,
       const GF3D2<add_const_t<T> > &gf_Ayy_,
       const GF3D2<add_const_t<T> > &gf_Ayz_,
       const GF3D2<add_const_t<T> > &gf_Azz_, const vect<int, 3> &I)
      : mat3{gf_Axx_(I), gf_Axy_(I), gf_Axz_(I),
             gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat3(const GF3D2<remove_const_t<T> > &gf_Axx_,
       const GF3D2<remove_const_t<T> > &gf_Axy_,
       const GF3D2<remove_const_t<T> > &gf_Axz_,
       const GF3D2<remove_const_t<T> > &gf_Ayy_,
       const GF3D2<remove_const_t<T> > &gf_Ayz_,
       const GF3D2<remove_const_t<T> > &gf_Azz_, const vect<int, 3> &I)
      : mat3{gf_Axx_(I), gf_Axy_(I), gf_Axz_(I),
             gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  template <typename F, typename = result_of_t<F(int, int)> >
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3(const F &f)
      : elts{f(0, 0), f(0, 1), f(0, 2), f(1, 1), f(1, 2), f(2, 2)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void
  store(const GF3D2<T> &gf_Axx_, const GF3D2<T> &gf_Axy_,
        const GF3D2<T> &gf_Axz_, const GF3D2<T> &gf_Ayy_,
        const GF3D2<T> &gf_Ayz_, const GF3D2<T> &gf_Azz_,
        const vect<int, 3> &I) const {
    const auto &A = *this;
    gf_Axx_(I) = A(0, 0);
    gf_Axy_(I) = A(0, 1);
    gf_Axz_(I) = A(0, 2);
    gf_Ayy_(I) = A(1, 1);
    gf_Ayz_(I) = A(1, 2);
    gf_Azz_(I) = A(2, 2);
  }

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST const T &
  operator()(int i, int j) const {
    return elts[ind(i, j)];
  }
  //  T &operator()(int i, int j) { return
  // elts[symind(i, j)]; }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T &
  operator()(int i, int j) {
    return elts[ind(i, j)];
  }

  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      mat3<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >,
           dnup1, dnup2>
      operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I),
            elts[3](I), elts[4](I), elts[5](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      mat3<remove_cv_t<
               remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
           dnup1, dnup2>
      operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[3](layout, I), elts[4](layout, I), elts[5](layout, I)};
  }

  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3<T, dnup1, dnup2>
      operator+(const mat3<T, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3<T, dnup1, dnup2>
      operator-(const mat3<T, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3<T, dnup1, dnup2>
      operator+(const mat3<T, dnup1, dnup2> &x,
                const mat3<T, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3<T, dnup1, dnup2>
      operator-(const mat3<T, dnup1, dnup2> &x,
                const mat3<T, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3<T, dnup1, dnup2>
      operator*(const T &a, const mat3<T, dnup1, dnup2> &x) {
    return {a * x.elts};
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator==(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return equal_to<vect<T, 6> >()(x.elts, y.elts);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator!=(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return !(x == y);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  maxabs() const {
    return elts.maxabs();
  }

  friend struct norm1<mat3>;

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T det() const {
    const auto &A = *this;
    return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) -
           A(1, 0) * (A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) +
           A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
  }

  constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3<T, !dnup1, !dnup2>
      inv(const T detA) const {
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

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  trace(const mat3<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    return sum2([&](int x, int y) { return gu(x, y) * A(x, y); });
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3 trace_free(
      const mat3<T, dnup1, dnup2> &g, const mat3<T, !dnup1, !dnup2> &gu) const {
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
  constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3<T, dnup1, dnup2>
      operator()() const {
    return mat3<T, dnup1, dnup2>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup1, Weyl::dnup_t dnup2>
struct zero<Weyl::mat3<T, dnup1, dnup2> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      Weyl::mat3<T, dnup1, dnup2>
      operator()() const {
    return Weyl::mat3<T, dnup1, dnup2>(zero<vect<T, 6> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct norm1<mat3<T, dnup1, dnup2> > {
  typedef typename norm1<vect<T, 6> >::result_type result_type;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST result_type
  operator()(const mat3<T, dnup1, dnup2> &x) const {
    return norm1<vect<T, 6> >()(x.elts);
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup3>
constexpr
    CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat3<T, dnup1, dnup2>
    mul(const mat3<T, dnup1, dnup3> &A, const mat3<T, !dnup3, dnup2> &B) {
  // C[a,b] = A[a,c] B[c,b]
  return mat3<T, dnup1, dnup2>([&](int a, int b) {
    return sum1([&](int x) { return A(a, x) * B(x, b); });
  });
}

template <typename F,
          typename R = remove_cv_t<remove_reference_t<result_of_t<F(int)> > > >
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
sum1(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 3; ++x)
    s += f(x);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int)> > > >
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
sum2(const F &f) {
  R s = zero<R>()();
  for (int x = 0; x < 3; ++x)
    for (int y = 0; y < 3; ++y)
      s += f(x, y);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int, int)> > > >
constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST R
sum3(const F &f) {
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

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 4);
#endif
    return n;
  }

public:
  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4()
      : elts{nan<vect<T, 4> >()()} {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec4(vect<T, 4> elts)
      : elts(move(elts)) {}

  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec4(T vt, T vx, T vy, T vz)
      : elts(make_tuple(move(vt), move(vx), move(vy), move(vz))) {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec4(initializer_list<T> v)
      : elts(v) {}
  // constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4(const
  // vector<T> &v) : elts(v) {} constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST vec4(vector<T> &&v) : elts(move(v)) {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec4(const GF3D2<add_const_t<T> > &gf_vt_,
       const GF3D2<add_const_t<T> > &gf_vx_,
       const GF3D2<add_const_t<T> > &gf_vy_,
       const GF3D2<add_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec4{gf_vt_(I), gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  vec4(const GF3D2<remove_const_t<T> > &gf_vt_,
       const GF3D2<remove_const_t<T> > &gf_vx_,
       const GF3D2<remove_const_t<T> > &gf_vy_,
       const GF3D2<remove_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec4{gf_vt_(I), gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  template <typename F, typename = result_of_t<F(int)> >
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4(const F &f)
      : elts{f(0), f(1), f(2), f(3)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void
  store(const GF3D2<T> &gf_vt_, const GF3D2<T> &gf_vx_, const GF3D2<T> &gf_vy_,
        const GF3D2<T> &gf_vz_, const vect<int, 3> &I) const {
    const auto &v = *this;
    gf_vt_(I) = v(0);
    gf_vx_(I) = v(1);
    gf_vy_(I) = v(2);
    gf_vz_(I) = v(3);
  }

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST const T &
  operator()(int i) const {
    return elts[ind(i)];
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T &
  operator()(int i) {
    return elts[ind(i)];
  }

  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<
      remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >, dnup>
  operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I), elts[3](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      vec4<remove_cv_t<
               remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 4>)> > >,
           dnup>
      operator()(const GF3D5layout &layout, const vect<int, 4> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[3](layout, I)};
  }

  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup>
      operator+(const vec4<T, dnup> &x) {
    return {+x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup>
      operator-(const vec4<T, dnup> &x) {
    return {-x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup>
      operator+(const vec4<T, dnup> &x, const vec4<T, dnup> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup>
      operator-(const vec4<T, dnup> &x, const vec4<T, dnup> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup>
      operator*(const T &a, const vec4<T, dnup> &x) {
    return {a * x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup>
      operator*(const vec4<T, dnup> &x, const T &a) {
    return {x.elts * a};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup>
      operator/(const vec4<T, dnup> &x, const T &a) {
    return {x.elts / a};
  }

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup> &
  operator*=(const T &a) {
    return *this = *this * a;
  }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup> &
  operator/=(const T &a) {
    return *this = *this / a;
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator==(const vec4<T, dnup> &x, const vec4<T, dnup> &y) {
    return equal_to<vect<T, 4> >()(x.elts, y.elts);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator!=(const vec4<T, dnup> &x, const vec4<T, dnup> &y) {
    return !(x == y);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<bool, dnup>
  isnan() const {
    using std::isnan;
    return vec4<bool, dnup>(isnan(elts));
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  any() const {
    return elts.any();
  }
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  maxabs() const {
    return elts.maxabs();
  }

  friend struct norm1<vec4>;

  friend ostream &operator<<(ostream &os, const vec4<T, dnup> &v) {
    return os << "[" << v(0) << "," << v(1) << "," << v(2) << "," << v(3)
              << "]";
  }
};

template <typename T, dnup_t dnup> struct nan<vec4<T, dnup> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST vec4<T, dnup>
  operator()() const {
    return vec4<T, dnup>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup> struct zero<Weyl::vec4<T, dnup> > {
  constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST Weyl::vec4<T, dnup>
      operator()() const {
    return Weyl::vec4<T, dnup>(zero<vect<T, 4> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup> struct norm1<vec4<T, dnup> > {
  typedef typename norm1<vect<T, 4> >::result_type result_type;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST result_type
  operator()(const vec4<T, dnup> &x) const {
    return norm1<vect<T, 4> >()(x.elts);
  }
};

////////////////////////////////////////////////////////////////////////////////

// Symmetric 4-matrix
template <typename T, dnup_t dnup1, dnup_t dnup2> class mat4 {
  static_assert(dnup1 == dnup2, "");

  vect<T, 10> elts;

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  symind(const int i, const int j) {
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
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  ind(const int i, const int j) {
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

  // nvcc doesn't handle these constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST expressions
#ifndef __CUDACC__
  static_assert(ind(1, 0) == ind(0, 1), "");
  static_assert(ind(2, 0) == ind(0, 2), "");
  static_assert(ind(3, 0) == ind(0, 3), "");
  static_assert(ind(2, 1) == ind(1, 2), "");
  static_assert(ind(3, 1) == ind(1, 3), "");
  static_assert(ind(3, 2) == ind(2, 3), "");
#endif

public:
  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4()
      : elts{nan<vect<T, 10> >()()} {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat4(vect<T, 10> elts)
      : elts(move(elts)) {}

  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat4(T Att, T Atx, T Aty, T Atz, T Axx, T Axy, T Axz, T Ayy, T Ayz, T Azz)
      : elts(make_tuple(move(Att), move(Atx), move(Aty), move(Atz), move(Axx),
                        move(Axy), move(Axz), move(Ayy), move(Ayz),
                        move(Azz))) {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  mat4(initializer_list<T> A)
      : elts(A) {}
  // constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4(const
  // vector<T> &A) : elts(A) {} constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST mat4(vector<T> &&A) : elts(move(A)) {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
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

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
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
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4(const F &f)
      : elts{f(0, 0), f(0, 1), f(0, 2), f(0, 3), f(1, 1),
             f(1, 2), f(1, 3), f(2, 2), f(2, 3), f(3, 3)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void store(
      const GF3D2<T> &gf_Att_, const GF3D2<T> &gf_Atx_, const GF3D2<T> &gf_Aty_,
      const GF3D2<T> &gf_Atz_, const GF3D2<T> &gf_Axx_, const GF3D2<T> &gf_Axy_,
      const GF3D2<T> &gf_Axz_, const GF3D2<T> &gf_Ayy_, const GF3D2<T> &gf_Ayz_,
      const GF3D2<T> &gf_Azz_, const vect<int, 3> &I) const {
    const auto &A = *this;
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

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST const T &
  operator()(int i, int j) const {
    return elts[ind(i, j)];
  }
  //  T &operator()(int i, int j) { return
  // elts[symind(i, j)]; }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T &
  operator()(int i, int j) {
    return elts[ind(i, j)];
  }

  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      mat4<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 10>)> > >,
           dnup1, dnup2>
      operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I), elts[3](I), elts[4](I),
            elts[5](I), elts[6](I), elts[7](I), elts[8](I), elts[9](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      mat4<remove_cv_t<
               remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
           dnup1, dnup2>
      operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[3](layout, I), elts[4](layout, I), elts[5](layout, I),
            elts[6](layout, I), elts[7](layout, I), elts[8](layout, I),
            elts[9](layout, I)};
  }

  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4<T, dnup1, dnup2>
      operator+(const mat4<T, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4<T, dnup1, dnup2>
      operator-(const mat4<T, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4<T, dnup1, dnup2>
      operator+(const mat4<T, dnup1, dnup2> &x,
                const mat4<T, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4<T, dnup1, dnup2>
      operator-(const mat4<T, dnup1, dnup2> &x,
                const mat4<T, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4<T, dnup1, dnup2>
      operator*(const T &a, const mat4<T, dnup1, dnup2> &x) {
    return {a * x.elts};
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator==(const mat4<T, dnup1, dnup2> &x, const mat4<T, dnup1, dnup2> &y) {
    return equal_to<vect<T, 10> >()(x.elts, y.elts);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator!=(const mat4<T, dnup1, dnup2> &x, const mat4<T, dnup1, dnup2> &y) {
    return !(x == y);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  maxabs() const {
    return elts.maxabs();
  }

  friend struct norm1<mat4>;

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T det() const {
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

  constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4<T, !dnup1, !dnup2>
      inv(const T detA) const {
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

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  trace(const mat4<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    return sum2([&](int x, int y) { return gu(x, y) * A(x, y); });
  }

#if 0
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST  mat4 trace_free(
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
  constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4<T, dnup1, dnup2>
      operator()() const {
    return mat4<T, dnup1, dnup2>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup1, Weyl::dnup_t dnup2>
struct zero<Weyl::mat4<T, dnup1, dnup2> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      Weyl::mat4<T, dnup1, dnup2>
      operator()() const {
    return Weyl::mat4<T, dnup1, dnup2>(zero<vect<T, 10> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct norm1<mat4<T, dnup1, dnup2> > {
  typedef typename norm1<vect<T, 10> >::result_type result_type;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST result_type
  operator()(const mat4<T, dnup1, dnup2> &x) const {
    return norm1<vect<T, 10> >()(x.elts);
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup4>
constexpr
    CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST mat4<T, dnup1, dnup2>
    mul(const mat4<T, dnup1, dnup4> &A, const mat4<T, !dnup4, dnup2> &B) {
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

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  asymind(const int i, const int j) {
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
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  ind(const int i, const int j) {
    return asymind(min(i, j), max(i, j));
  }
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  sign(const int i, const int j) {
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

  // nvcc doesn't handle these constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST expressions
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
  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4()
      : elts{nan<vect<T, 6> >()()} {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  amat4(vect<T, 6> elts)
      : elts(move(elts)) {}

  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  amat4(T Atx, T Aty, T Atz, T Axy, T Axz, T Ayz)
      : elts(make_tuple(move(Atx), move(Aty), move(Atz), move(Axy), move(Axz),
                        move(Ayz))) {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  amat4(initializer_list<T> A)
      : elts(A) {}
  // constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  // amat4(const vector<T> &A) : elts(A) {} constexpr
  // CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST amat4(vector<T> &&A) : elts(move(A)) {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  amat4(const GF3D2<add_const_t<T> > &gf_Atx_,
        const GF3D2<add_const_t<T> > &gf_Aty_,
        const GF3D2<add_const_t<T> > &gf_Atz_,
        const GF3D2<add_const_t<T> > &gf_Axy_,
        const GF3D2<add_const_t<T> > &gf_Axz_,
        const GF3D2<add_const_t<T> > &gf_Ayz_, const vect<int, 3> &I)
      : amat4{gf_Atx_(I), gf_Aty_(I), gf_Atz_(I),
              gf_Axy_(I), gf_Axz_(I), gf_Ayz_(I)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  amat4(const GF3D2<remove_const_t<T> > &gf_Atx_,
        const GF3D2<remove_const_t<T> > &gf_Aty_,
        const GF3D2<remove_const_t<T> > &gf_Atz_,
        const GF3D2<remove_const_t<T> > &gf_Axy_,
        const GF3D2<remove_const_t<T> > &gf_Axz_,
        const GF3D2<remove_const_t<T> > &gf_Ayz_, const vect<int, 3> &I)
      : amat4{gf_Atx_(I), gf_Aty_(I), gf_Atz_(I),
              gf_Axy_(I), gf_Axz_(I), gf_Ayz_(I)} {}

  template <typename F, typename = result_of_t<F(int, int)> >
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4(const F &f)
      : elts{f(0, 1), f(0, 2), f(0, 3), f(1, 2), f(1, 3), f(2, 3)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void
  store(const GF3D2<T> &gf_Atx_, const GF3D2<T> &gf_Aty_,
        const GF3D2<T> &gf_Atz_, const GF3D2<T> &gf_Axy_,
        const GF3D2<T> &gf_Axz_, const GF3D2<T> &gf_Ayz_,
        const vect<int, 3> &I) const {
    const auto &A = *this;
    gf_Atx_(I) = A(0, 1);
    gf_Aty_(I) = A(0, 2);
    gf_Atz_(I) = A(0, 3);
    gf_Axy_(I) = A(1, 2);
    gf_Axz_(I) = A(1, 3);
    gf_Ayz_(I) = A(2, 3);
  }

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST const T
  operator()(int i, int j) const {
    const int s = sign(i, j);
    if (s == 0)
      return zero<T>()();
    return s * elts[ind(i, j)];
  }
  //  T &operator()(int i, int j) { return
  // elts[symind(i, j)]; }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T &
  operator()(int i, int j) {
    const int s = sign(i, j);
    assert(s == 1);
    return elts[ind(i, j)];
  }

  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      amat4<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 6>)> > >,
            dnup1, dnup2>
      operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I),
            elts[3](I), elts[4](I), elts[5](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      amat4<remove_cv_t<remove_reference_t<
                result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
            dnup1, dnup2>
      operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[3](layout, I), elts[4](layout, I), elts[5](layout, I)};
  }

  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4<T, dnup1, dnup2>
      operator+(const amat4<T, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4<T, dnup1, dnup2>
      operator-(const amat4<T, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4<T, dnup1, dnup2>
      operator+(const amat4<T, dnup1, dnup2> &x,
                const amat4<T, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4<T, dnup1, dnup2>
      operator-(const amat4<T, dnup1, dnup2> &x,
                const amat4<T, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4<T, dnup1, dnup2>
      operator*(const T &a, const amat4<T, dnup1, dnup2> &x) {
    return {a * x.elts};
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator==(const amat4<T, dnup1, dnup2> &x, const amat4<T, dnup1, dnup2> &y) {
    return equal_to<vect<T, 6> >()(x.elts, y.elts);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator!=(const amat4<T, dnup1, dnup2> &x, const amat4<T, dnup1, dnup2> &y) {
    return !(x == y);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  maxabs() const {
    return elts.maxabs();
  }

  friend struct norm1<amat4>;

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T det() const {
    const auto &A = *this;
    return pow(A(0, 3), 2) * pow(A(1, 2), 2) -
           2 * A(0, 2) * A(0, 3) * A(1, 2) * A(1, 3) +
           pow(A(0, 2), 2) * pow(A(1, 3), 2) +
           2 * A(0, 1) * A(0, 3) * A(1, 2) * A(2, 3) -
           2 * A(0, 1) * A(0, 2) * A(1, 3) * A(2, 3) +
           pow(A(0, 1), 2) * pow(A(2, 3), 2);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      amat4<T, !dnup1, !dnup2>
      inv(const T detA) const {
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

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  trace(const amat4<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    return sum2([&](int x, int y) { return gu(x, y) * A(x, y); });
  }

#if 0
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4
  trace_free(const amat4<T, dnup1, dnup2> &g,
             const amat4<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    const T trA = A.trace(gu);
    return amat4([&](int a, int b) { return A(a, b) - trA / 2 * g(a, b); });
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
  constexpr
      CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4<T, dnup1, dnup2>
      operator()() const {
    return amat4<T, dnup1, dnup2>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup1, Weyl::dnup_t dnup2>
struct zero<Weyl::amat4<T, dnup1, dnup2> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      Weyl::amat4<T, dnup1, dnup2>
      operator()() const {
    return Weyl::amat4<T, dnup1, dnup2>(zero<vect<T, 6> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct norm1<amat4<T, dnup1, dnup2> > {
  typedef typename norm1<vect<T, 6> >::result_type result_type;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST result_type
  operator()(const amat4<T, dnup1, dnup2> &x) const {
    return norm1<vect<T, 6> >()(x.elts);
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup4>
constexpr
    CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST amat4<T, dnup1, dnup2>
    mul(const amat4<T, dnup1, dnup4> &A, const amat4<T, !dnup4, dnup2> &B) {
  // C[a,b] = A[a,c] B[c,b]
  return amat4<T, dnup1, dnup2>([&](int a, int b) {
    return sum1([&](int x) { return A(a, x) * B(x, b); });
  });
}

////////////////////////////////////////////////////////////////////////////////

// Riemann-tensor like object (4 dimensiona, rank 4)
template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup3, dnup_t dnup4>
class rten4 {
  static_assert(dnup1 == dnup2 && dnup1 == dnup3 && dnup1 == dnup4, "");

  // We omit the first Bianchi identity
  // R_abcd = - R_abdc
  // R_abcd = - R_bacd
  // R_abcd = + R_cdab
  vect<T, 21> elts;

  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  rsymind(const int i, const int j, const int k, const int l) {
#ifdef CCTK_DEBUG
    assert(i >= 0 && i < j && j < 4);
    assert(k >= 0 && k < l && l < 4);
#endif
    const int ij = i * (5 - i) / 2 + j - 1; // asymind D=4
    const int kl = k * (5 - k) / 2 + l - 1; // asymind D=4
#ifdef CCTK_DEBUG
    assert(ij >= 0 && ij <= kl && kl < 6);
#endif
    const int n = ij * (11 - ij) / 2 + kl; // symind D=6

    // i j   k l   ij   kl    n
    // 0 1   0 1    0    0    0
    // 0 1   0 2    0    1    1
    // 0 1   0 3    0    2    2
    // 0 1   1 2    0    3    3
    // 0 1   1 3    0    4    4
    // 0 1   2 3    0    5    5
    // 0 2   0 2    1    1    6
    // 0 2   0 3    1    2    7
    // 0 2   1 2    1    3    8
    // 0 2   1 3    1    4    9
    // 0 2   2 3    1    5   10
    // 0 3   0 3    2    2   11
    // 0 3   1 2    2    3   12
    // 0 3   1 3    2    4   13
    // 0 3   2 3    2    5   14
    // 1 2   1 2    3    3   15
    // 1 2   1 3    3    4   16
    // 1 2   2 3    3    5   17
    // 1 3   1 3    4    4   18
    // 1 3   2 3    4    5   19
    // 2 3   2 3    5    5   20
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 21);
#endif
    return n;
  }
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  ind(const int i, const int j, const int k, const int l) {
    const int i1 = min(i, j);
    const int j1 = max(i, j);
    const int k1 = min(k, l);
    const int l1 = max(k, l);
    const bool noswap = 4 * i1 + j1 <= 4 * k1 + l1;
    const int i2 = noswap ? i1 : k1;
    const int j2 = noswap ? j1 : l1;
    const int k2 = noswap ? k1 : i1;
    const int l2 = noswap ? l1 : j1;
    return rsymind(i2, j2, k2, l2);
  }
  static constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int
  sign(const int i, const int j, const int k, const int l) {
    if (i == j || k == l)
      return 0;
    return (i < j) ^ (k < l) ? 1 : -1;
  }

  static_assert(rsymind(0, 1, 0, 1) == 0, "");
  static_assert(rsymind(0, 1, 0, 2) == 1, "");
  static_assert(rsymind(0, 1, 0, 3) == 2, "");
  static_assert(rsymind(0, 1, 1, 2) == 3, "");
  static_assert(rsymind(0, 1, 1, 3) == 4, "");
  static_assert(rsymind(0, 1, 2, 3) == 5, "");
  static_assert(rsymind(0, 2, 0, 2) == 6, "");
  static_assert(rsymind(0, 2, 0, 3) == 7, "");
  static_assert(rsymind(0, 2, 1, 2) == 8, "");
  static_assert(rsymind(0, 2, 1, 3) == 9, "");
  static_assert(rsymind(0, 2, 2, 3) == 10, "");
  static_assert(rsymind(0, 3, 0, 3) == 11, "");
  static_assert(rsymind(0, 3, 1, 2) == 12, "");
  static_assert(rsymind(0, 3, 1, 3) == 13, "");
  static_assert(rsymind(0, 3, 2, 3) == 14, "");
  static_assert(rsymind(1, 2, 1, 2) == 15, "");
  static_assert(rsymind(1, 2, 1, 3) == 16, "");
  static_assert(rsymind(1, 2, 2, 3) == 17, "");
  static_assert(rsymind(1, 3, 1, 3) == 18, "");
  static_assert(rsymind(1, 3, 2, 3) == 19, "");
  static_assert(rsymind(2, 3, 2, 3) == 20, "");

  // nvcc doesn't handle these constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST expressions
#ifndef __CUDACC__
  static_assert(ind(0, 1, 1, 0) == ind(0, 1, 0, 1), "");
  static_assert(ind(1, 0, 0, 1) == ind(0, 1, 0, 1), "");
  static_assert(ind(1, 0, 1, 0) == ind(0, 1, 0, 1), "");
  // .. there are too many to write out
#endif

  static_assert(sign(0, 1, 1, 0) == -sign(0, 1, 0, 1), "");
  static_assert(sign(1, 0, 0, 1) == -sign(0, 1, 0, 1), "");
  static_assert(sign(1, 0, 1, 0) == +sign(0, 1, 0, 1), "");
  // ... there are too many to write out

  static_assert(sign(0, 0, 0, 0) == 0, "");
  static_assert(sign(0, 1, 0, 0) == 0, "");
  static_assert(sign(0, 0, 0, 1) == 0, "");
  // ... there are too many to write out

public:
  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST rten4()
      : elts{nan<vect<T, 21> >()()} {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  rten4(vect<T, 21> elts)
      : elts(move(elts)) {}

  explicit constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  rten4(T Rm0101, T Rm0102, T Rm0103, T Rm0112, T Rm0113, T Rm0123, T Rm0202,
        T Rm0203, T Rm0212, T Rm0213, T Rm0223, T Rm0303, T Rm0312, T Rm0313,
        T Rm0323, T Rm1212, T Rm1213, T Rm1223, T Rm1313, T Rm1323, T Rm2323)
      : elts(make_tuple(move(Rm0101), move(Rm0102), move(Rm0103), move(Rm0112),
                        move(Rm0113), move(Rm0123), move(Rm0202), move(Rm0203),
                        move(Rm0212), move(Rm0213), move(Rm0223), move(Rm0303),
                        move(Rm0312), move(Rm0313), move(Rm0323), move(Rm1212),
                        move(Rm1213), move(Rm1223), move(Rm1313), move(Rm1323),
                        move(Rm2323))) {}

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  rten4(initializer_list<T> A)
      : elts(A) {}
  // constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  // rten4(const vector<T> &A) : elts(A) {} constexpr
  // CCTK_ATTRIBUTE_ALWAYS_INLINE
  // CCTK_DEVICE CCTK_HOST rten4(vector<T> &&A) : elts(move(A)) {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  rten4(const GF3D2<add_const_t<T> > &gf_Rm0101_,
        const GF3D2<add_const_t<T> > &gf_Rm0102_,
        const GF3D2<add_const_t<T> > &gf_Rm0103_,
        const GF3D2<add_const_t<T> > &gf_Rm0112_,
        const GF3D2<add_const_t<T> > &gf_Rm0113_,
        const GF3D2<add_const_t<T> > &gf_Rm0123_,
        const GF3D2<add_const_t<T> > &gf_Rm0202_,
        const GF3D2<add_const_t<T> > &gf_Rm0203_,
        const GF3D2<add_const_t<T> > &gf_Rm0212_,
        const GF3D2<add_const_t<T> > &gf_Rm0213_,
        const GF3D2<add_const_t<T> > &gf_Rm0223_,
        const GF3D2<add_const_t<T> > &gf_Rm0303_,
        const GF3D2<add_const_t<T> > &gf_Rm0312_,
        const GF3D2<add_const_t<T> > &gf_Rm0313_,
        const GF3D2<add_const_t<T> > &gf_Rm0323_,
        const GF3D2<add_const_t<T> > &gf_Rm1212_,
        const GF3D2<add_const_t<T> > &gf_Rm1213_,
        const GF3D2<add_const_t<T> > &gf_Rm1223_,
        const GF3D2<add_const_t<T> > &gf_Rm1313_,
        const GF3D2<add_const_t<T> > &gf_Rm1323_,
        const GF3D2<add_const_t<T> > &gf_Rm2323_, const vect<int, 3> &I)
      : rten4{gf_Rm0101_(I), gf_Rm0102_(I), gf_Rm0103_(I), gf_Rm0112_(I),
              gf_Rm0113_(I), gf_Rm0123_(I), gf_Rm0202_(I), gf_Rm0203_(I),
              gf_Rm0212_(I), gf_Rm0213_(I), gf_Rm0223_(I), gf_Rm0303_(I),
              gf_Rm0312_(I), gf_Rm0313_(I), gf_Rm0323_(I), gf_Rm1212_(I),
              gf_Rm1213_(I), gf_Rm1223_(I), gf_Rm1313_(I), gf_Rm1323_(I),
              gf_Rm2323_(I)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
  rten4(const GF3D2<remove_const_t<T> > &gf_Rm0101_,
        const GF3D2<remove_const_t<T> > &gf_Rm0102_,
        const GF3D2<remove_const_t<T> > &gf_Rm0103_,
        const GF3D2<remove_const_t<T> > &gf_Rm0112_,
        const GF3D2<remove_const_t<T> > &gf_Rm0113_,
        const GF3D2<remove_const_t<T> > &gf_Rm0123_,
        const GF3D2<remove_const_t<T> > &gf_Rm0202_,
        const GF3D2<remove_const_t<T> > &gf_Rm0203_,
        const GF3D2<remove_const_t<T> > &gf_Rm0212_,
        const GF3D2<remove_const_t<T> > &gf_Rm0213_,
        const GF3D2<remove_const_t<T> > &gf_Rm0223_,
        const GF3D2<remove_const_t<T> > &gf_Rm0303_,
        const GF3D2<remove_const_t<T> > &gf_Rm0312_,
        const GF3D2<remove_const_t<T> > &gf_Rm0313_,
        const GF3D2<remove_const_t<T> > &gf_Rm0323_,
        const GF3D2<remove_const_t<T> > &gf_Rm1212_,
        const GF3D2<remove_const_t<T> > &gf_Rm1213_,
        const GF3D2<remove_const_t<T> > &gf_Rm1223_,
        const GF3D2<remove_const_t<T> > &gf_Rm1313_,
        const GF3D2<remove_const_t<T> > &gf_Rm1323_,
        const GF3D2<remove_const_t<T> > &gf_Rm2323_, const vect<int, 3> &I)
      : rten4{gf_Rm0101_(I), gf_Rm0102_(I), gf_Rm0103_(I), gf_Rm0112_(I),
              gf_Rm0113_(I), gf_Rm0123_(I), gf_Rm0202_(I), gf_Rm0203_(I),
              gf_Rm0212_(I), gf_Rm0213_(I), gf_Rm0223_(I), gf_Rm0303_(I),
              gf_Rm0312_(I), gf_Rm0313_(I), gf_Rm0323_(I), gf_Rm1212_(I),
              gf_Rm1213_(I), gf_Rm1223_(I), gf_Rm1313_(I), gf_Rm1323_(I),
              gf_Rm2323_(I)} {}

  template <typename F, typename = result_of_t<F(int, int, int, int)> >
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST rten4(const F &f)
      : elts{f(0, 1, 0, 1), f(0, 1, 0, 2), f(0, 1, 0, 3), f(0, 1, 1, 2),
             f(0, 1, 1, 3), f(0, 1, 2, 3), f(0, 2, 0, 2), f(0, 2, 0, 3),
             f(0, 2, 1, 2), f(0, 2, 1, 3), f(0, 2, 2, 3), f(0, 3, 0, 3),
             f(0, 3, 1, 2), f(0, 3, 1, 3), f(0, 3, 2, 3), f(1, 2, 1, 2),
             f(1, 2, 1, 3), f(1, 2, 2, 3), f(1, 3, 1, 3), f(1, 3, 2, 3),
             f(2, 3, 2, 3)} {}

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST void
  store(const GF3D2<T> &gf_Rm0101_, const GF3D2<T> &gf_Rm0102_,
        const GF3D2<T> &gf_Rm0103_, const GF3D2<T> &gf_Rm0112_,
        const GF3D2<T> &gf_Rm0113_, const GF3D2<T> &gf_Rm0123_,
        const GF3D2<T> &gf_Rm0202_, const GF3D2<T> &gf_Rm0203_,
        const GF3D2<T> &gf_Rm0212_, const GF3D2<T> &gf_Rm0213_,
        const GF3D2<T> &gf_Rm0223_, const GF3D2<T> &gf_Rm0303_,
        const GF3D2<T> &gf_Rm0312_, const GF3D2<T> &gf_Rm0313_,
        const GF3D2<T> &gf_Rm0323_, const GF3D2<T> &gf_Rm1212_,
        const GF3D2<T> &gf_Rm1213_, const GF3D2<T> &gf_Rm1223_,
        const GF3D2<T> &gf_Rm1313_, const GF3D2<T> &gf_Rm1323_,
        const GF3D2<T> &gf_Rm2323_, const vect<int, 3> &I) const {
    const auto &Rm = *this;
    gf_Rm0101_(I) = Rm(0, 1, 0, 1);
    gf_Rm0102_(I) = Rm(0, 1, 0, 2);
    gf_Rm0103_(I) = Rm(0, 1, 0, 3);
    gf_Rm0112_(I) = Rm(0, 1, 1, 2);
    gf_Rm0113_(I) = Rm(0, 1, 1, 3);
    gf_Rm0123_(I) = Rm(0, 1, 2, 3);
    gf_Rm0202_(I) = Rm(0, 2, 0, 2);
    gf_Rm0203_(I) = Rm(0, 2, 0, 3);
    gf_Rm0212_(I) = Rm(0, 2, 1, 2);
    gf_Rm0213_(I) = Rm(0, 2, 1, 3);
    gf_Rm0223_(I) = Rm(0, 2, 2, 3);
    gf_Rm0303_(I) = Rm(0, 3, 0, 3);
    gf_Rm0312_(I) = Rm(0, 3, 1, 2);
    gf_Rm0313_(I) = Rm(0, 3, 1, 3);
    gf_Rm0323_(I) = Rm(0, 3, 2, 3);
    gf_Rm1212_(I) = Rm(1, 2, 1, 2);
    gf_Rm1213_(I) = Rm(1, 2, 1, 3);
    gf_Rm1223_(I) = Rm(1, 2, 2, 3);
    gf_Rm1313_(I) = Rm(1, 3, 1, 3);
    gf_Rm1323_(I) = Rm(1, 3, 2, 3);
    gf_Rm2323_(I) = Rm(2, 3, 2, 3);
  }

  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST const T
  operator()(int i, int j, int k, int l) const {
    const int s = sign(i, j, k, l);
    if (s == 0)
      return zero<T>()();
    return s * elts[ind(i, j, k, l)];
  }
  //  T &operator()(int i, int j, int k, int l) { return
  // elts[symind(i, j, k, l)]; }
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T &
  operator()(int i, int j, int k, int l) {
    const int s = sign(i, j, k, l);
    assert(s == 1);
    return elts[ind(i, j, k, l)];
  }

  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      rten4<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 21>)> > >,
            dnup1, dnup2, dnup3, dnup4>
      operator()(const vect<int, 3> &I) const {
    return {elts[0](I),  elts[1](I),  elts[2](I),  elts[3](I),  elts[4](I),
            elts[5](I),  elts[6](I),  elts[7](I),  elts[8](I),  elts[9](I),
            elts[10](I), elts[11](I), elts[12](I), elts[13](I), elts[14](I),
            elts[15](I), elts[16](I), elts[17](I), elts[18](I), elts[19](I),
            elts[20](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      rten4<remove_cv_t<remove_reference_t<
                result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
            dnup1, dnup2, dnup3, dnup4>
      operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I),  elts[1](layout, I),  elts[2](layout, I),
            elts[3](layout, I),  elts[4](layout, I),  elts[5](layout, I),
            elts[6](layout, I),  elts[7](layout, I),  elts[8](layout, I),
            elts[9](layout, I),  elts[10](layout, I), elts[11](layout, I),
            elts[12](layout, I), elts[13](layout, I), elts[14](layout, I),
            elts[15](layout, I), elts[16](layout, I), elts[17](layout, I),
            elts[18](layout, I), elts[19](layout, I), elts[20](layout, I)};
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      rten4<T, dnup1, dnup2, dnup3, dnup4>
      operator+(const rten4<T, dnup1, dnup2, dnup3, dnup4> &x) {
    return {+x.elts};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      rten4<T, dnup1, dnup2, dnup3, dnup4>
      operator-(const rten4<T, dnup1, dnup2, dnup3, dnup4> &x) {
    return {-x.elts};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      rten4<T, dnup1, dnup2, dnup3, dnup4>
      operator+(const rten4<T, dnup1, dnup2, dnup3, dnup4> &x,
                const rten4<T, dnup1, dnup2, dnup3, dnup4> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      rten4<T, dnup1, dnup2, dnup3, dnup4>
      operator-(const rten4<T, dnup1, dnup2, dnup3, dnup4> &x,
                const rten4<T, dnup1, dnup2, dnup3, dnup4> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      rten4<T, dnup1, dnup2, dnup3, dnup4>
      operator*(const T &a, const rten4<T, dnup1, dnup2, dnup3, dnup4> &x) {
    return {a * x.elts};
  }

  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator==(const rten4<T, dnup1, dnup2, dnup3, dnup4> &x,
             const rten4<T, dnup1, dnup2, dnup3, dnup4> &y) {
    return equal_to<vect<T, 21> >()(x.elts, y.elts);
  }
  friend constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST bool
  operator!=(const rten4<T, dnup1, dnup2, dnup3, dnup4> &x,
             const rten4<T, dnup1, dnup2, dnup3, dnup4> &y) {
    return !(x == y);
  }

  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
  maxabs() const {
    return elts.maxabs();
  }

  friend struct norm1<rten4>;

  friend ostream &operator<<(ostream &os,
                             const rten4<T, dnup1, dnup2, dnup3, dnup4> &Rm) {
    return os << "[[" << Rm(0, 1, 0, 1) << "," << Rm(0, 1, 0, 2) << ","
              << Rm(0, 1, 0, 3) << "," << Rm(0, 1, 1, 2) << ","
              << Rm(0, 1, 1, 3) << "," << Rm(0, 1, 2, 3) << "],["
              << Rm(0, 2, 0, 2) << "," << Rm(0, 2, 0, 3) << ","
              << Rm(0, 2, 1, 2) << "," << Rm(0, 2, 1, 3) << ","
              << Rm(0, 2, 2, 3) << "],[" << Rm(0, 3, 0, 3) << ","
              << Rm(0, 3, 1, 2) << "," << Rm(0, 3, 1, 3) << ","
              << Rm(0, 3, 2, 3) << "],[" << Rm(1, 2, 1, 2) << ","
              << Rm(1, 2, 1, 3) << "," << Rm(1, 2, 2, 3) << "],["
              << Rm(1, 3, 1, 3) << "," << Rm(1, 3, 2, 3) << "],["
              << Rm(2, 3, 2, 3) << "]]";
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup3, dnup_t dnup4>
struct nan<rten4<T, dnup1, dnup2, dnup3, dnup4> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      rten4<T, dnup1, dnup2, dnup3, dnup4>
      operator()() const {
    return rten4<T, dnup1, dnup2, dnup3, dnup4>();
  }
};

} // namespace Weyl
namespace Arith {
template <typename T, Weyl::dnup_t dnup1, Weyl::dnup_t dnup2,
          Weyl::dnup_t dnup3, Weyl::dnup_t dnup4>
struct zero<Weyl::rten4<T, dnup1, dnup2, dnup3, dnup4> > {
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST
      Weyl::rten4<T, dnup1, dnup2, dnup3, dnup4>
      operator()() const {
    return Weyl::rten4<T, dnup1, dnup2, dnup3, dnup4>(zero<vect<T, 21> >()());
  }
};
} // namespace Arith
namespace Weyl {

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup3, dnup_t dnup4>
struct norm1<rten4<T, dnup1, dnup2, dnup3, dnup4> > {
  typedef typename norm1<vect<T, 21> >::result_type result_type;
  constexpr CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST result_type
  operator()(const rten4<T, dnup1, dnup2, dnup3, dnup4> &x) const {
    return norm1<vect<T, 21> >()(x.elts);
  }
};

} // namespace Weyl

#endif // #ifndef TENSOR_HXX
