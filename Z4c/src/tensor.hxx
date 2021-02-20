#ifndef TENSOR_HXX
#define TENSOR_HXX

#include <loop.hxx>

#include <cctk.h>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <sstream>
#include <tuple>
#include <type_traits>

#ifdef CCTK_DEBUG
#define Z4C_INLINE
#else
#define Z4C_INLINE CCTK_ATTRIBUTE_ALWAYS_INLINE
#endif

namespace Z4c {
using namespace Loop;
using namespace std;

template <typename T> Z4C_INLINE constexpr T pow2(const T x) { return x * x; }
template <typename T> Z4C_INLINE constexpr T pow3(const T x) {
  const T x2 = x * x;
  return x2 * x;
}
template <typename T> Z4C_INLINE constexpr T pow4(const T x) {
  const T x2 = x * x;
  return x2 * x2;
}
template <typename T> Z4C_INLINE constexpr T pow5(const T x) {
  const T x2 = x * x;
  const T x4 = x2 * x2;
  return x4 * x;
}
template <typename T> Z4C_INLINE constexpr T pow6(const T x) {
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

template <typename F, typename T>
constexpr Z4C_INLINE T fold(const F &f, const T &x) {
  return x;
}
template <typename F, typename T, typename... Ts>
constexpr Z4C_INLINE T fold(const F &f, const T &x0, const T &x1,
                            const Ts &...xs) {
  return fold(f, fold(f, x0, x1), xs...);
}

template <typename T> constexpr Z4C_INLINE T add() { return T(0); }
// template <typename T>
//  constexpr Z4C_INLINE T add(const T &x) {
//   return x;
// }
template <typename T, typename... Ts>
constexpr Z4C_INLINE T add(const T &x, const Ts &...xs) {
  return x + add(xs...);
}

template <typename F,
          typename R = remove_cv_t<remove_reference_t<result_of_t<F(int)> > > >
constexpr Z4C_INLINE R sum1(const F &f) {
  R s{0};
  for (int x = 0; x < 3; ++x)
    s += f(x);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int)> > > >
constexpr Z4C_INLINE R sum2(const F &f) {
  R s{0};
  for (int x = 0; x < 3; ++x)
    for (int y = 0; y < 3; ++y)
      s += f(x, y);
  return s;
}

template <typename F, typename R = remove_cv_t<
                          remove_reference_t<result_of_t<F(int, int, int)> > > >
constexpr Z4C_INLINE R sum3(const F &f) {
  R s{0};
  for (int x = 0; x < 3; ++x)
    for (int y = 0; y < 3; ++y)
      for (int z = 0; z < 3; ++z)
        s += f(x, y, z);
  return s;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct nan {
  constexpr Z4C_INLINE T operator()() const { return NAN; }
};
template <typename T, int D> struct nan<vect<T, D> > {
  constexpr Z4C_INLINE vect<T, D> operator()() const {
    return vect<T, D>::pure(nan<T>()());
  }
};

template <typename T> struct norm1 {
  typedef T result_type;
  constexpr Z4C_INLINE result_type operator()(const T &x) const {
    return abs(x);
  }
};
template <typename T, int D> struct norm1<vect<T, D> > {
  typedef typename norm1<T>::result_type result_type;
  constexpr Z4C_INLINE result_type operator()(const vect<T, D> &xs) const {
    typename norm1<T>::result_type r{0};
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

  static constexpr Z4C_INLINE int ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 3);
#endif
    return n;
  }

public:
  explicit constexpr Z4C_INLINE vec3() : elts{nan<vect<T, 3> >()()} {}

  constexpr Z4C_INLINE vec3(vect<T, 3> elts) : elts(move(elts)) {}

  explicit constexpr Z4C_INLINE vec3(T vx, T vy, T vz)
      : elts(make_tuple(move(vx), move(vy), move(vz))) {}

  constexpr Z4C_INLINE vec3(initializer_list<T> v) : elts(v) {}
  // constexpr Z4C_INLINE vec3(const vector<T> &v) : elts(v) {}
  // constexpr Z4C_INLINE vec3(vector<T> &&v) : elts(move(v)) {}

  Z4C_INLINE
  vec3(const GF3D2<add_const_t<T> > &gf_vx_,
       const GF3D2<add_const_t<T> > &gf_vy_,
       const GF3D2<add_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec3{gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}
  Z4C_INLINE
  vec3(const GF3D2<remove_const_t<T> > &gf_vx_,
       const GF3D2<remove_const_t<T> > &gf_vy_,
       const GF3D2<remove_const_t<T> > &gf_vz_, const vect<int, 3> &I)
      : vec3{gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  template <typename F, typename = result_of_t<F(int)> >
  constexpr Z4C_INLINE vec3(const F &f) : elts{f(0), f(1), f(2)} {}

  Z4C_INLINE void store(const GF3D2<T> &gf_vx_, const GF3D2<T> &gf_vy_,
                        const GF3D2<T> &gf_vz_, const vect<int, 3> &I) const {
    const auto &v = *this;
#ifdef CCTK_DEBUG
    assert(CCTK_isfinite(v(0)));
    assert(CCTK_isfinite(v(1)));
    assert(CCTK_isfinite(v(2)));
#endif
    gf_vx_(I) = v(0);
    gf_vy_(I) = v(1);
    gf_vz_(I) = v(2);
  }

  Z4C_INLINE const T &operator()(int i) const { return elts[ind(i)]; }
  Z4C_INLINE T &operator()(int i) { return elts[ind(i)]; }

  template <typename U = T>
  Z4C_INLINE
      vec3<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >,
           dnup>
      operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  Z4C_INLINE
      vec3<remove_cv_t<
               remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
           dnup>
      operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I)};
  }

  friend constexpr Z4C_INLINE vec3<T, dnup> operator+(const vec3<T, dnup> &x) {
    return {+x.elts};
  }
  friend constexpr Z4C_INLINE vec3<T, dnup> operator-(const vec3<T, dnup> &x) {
    return {-x.elts};
  }
  friend constexpr Z4C_INLINE vec3<T, dnup> operator+(const vec3<T, dnup> &x,
                                                      const vec3<T, dnup> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr Z4C_INLINE vec3<T, dnup> operator-(const vec3<T, dnup> &x,
                                                      const vec3<T, dnup> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr Z4C_INLINE vec3<T, dnup> operator*(const T &a,
                                                      const vec3<T, dnup> &x) {
    return {a * x.elts};
  }

  friend constexpr Z4C_INLINE bool operator==(const vec3<T, dnup> &x,
                                              const vec3<T, dnup> &y) {
    return equal_to<vect<T, 3> >()(x.elts, y.elts);
  }
  friend constexpr Z4C_INLINE bool operator!=(const vec3<T, dnup> &x,
                                              const vec3<T, dnup> &y) {
    return !(x == y);
  }

  constexpr Z4C_INLINE T maxabs() const { return elts.maxabs(); }

  friend struct norm1<vec3>;

  friend ostream &operator<<(ostream &os, const vec3<T, dnup> &v) {
    return os << "[" << v(0) << "," << v(1) << "," << v(2) << "]";
  }
};

template <typename T, dnup_t dnup> struct nan<vec3<T, dnup> > {
  constexpr Z4C_INLINE vec3<T, dnup> operator()() const {
    return vec3<T, dnup>();
  }
};

template <typename T, dnup_t dnup> struct norm1<vec3<T, dnup> > {
  typedef typename norm1<vect<T, 3> >::result_type result_type;
  constexpr Z4C_INLINE result_type operator()(const vec3<T, dnup> &x) const {
    return norm1<vect<T, 3> >()(x.elts);
  }
};

////////////////////////////////////////////////////////////////////////////////

// Symmetric 3-matrix
template <typename T, dnup_t dnup1, dnup_t dnup2> class mat3 {
  static_assert(dnup1 == dnup2, "");

  vect<T, 6> elts;

  static constexpr Z4C_INLINE int symind(const int i, const int j) {
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
    assert(n >= 0 && n < 6);
#endif
    return n;
  }
  static constexpr Z4C_INLINE int ind(const int i, const int j) {
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
  explicit constexpr Z4C_INLINE mat3() : elts{nan<vect<T, 6> >()()} {}

  constexpr Z4C_INLINE mat3(vect<T, 6> elts) : elts(move(elts)) {}

  explicit constexpr Z4C_INLINE mat3(T Axx, T Axy, T Axz, T Ayy, T Ayz, T Azz)
      : elts(make_tuple(move(Axx), move(Axy), move(Axz), move(Ayy), move(Ayz),
                        move(Azz))) {}

  constexpr Z4C_INLINE mat3(initializer_list<T> A) : elts(A) {}
  // constexpr Z4C_INLINE mat3(const vector<T> &A) : elts(A) {}
  // constexpr Z4C_INLINE mat3(vector<T> &&A) : elts(move(A)) {}

  Z4C_INLINE
  mat3(const GF3D2<add_const_t<T> > &gf_Axx_,
       const GF3D2<add_const_t<T> > &gf_Axy_,
       const GF3D2<add_const_t<T> > &gf_Axz_,
       const GF3D2<add_const_t<T> > &gf_Ayy_,
       const GF3D2<add_const_t<T> > &gf_Ayz_,
       const GF3D2<add_const_t<T> > &gf_Azz_, const vect<int, 3> &I)
      : mat3{gf_Axx_(I), gf_Axy_(I), gf_Axz_(I),
             gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}
  Z4C_INLINE
  mat3(const GF3D2<remove_const_t<T> > &gf_Axx_,
       const GF3D2<remove_const_t<T> > &gf_Axy_,
       const GF3D2<remove_const_t<T> > &gf_Axz_,
       const GF3D2<remove_const_t<T> > &gf_Ayy_,
       const GF3D2<remove_const_t<T> > &gf_Ayz_,
       const GF3D2<remove_const_t<T> > &gf_Azz_, const vect<int, 3> &I)
      : mat3{gf_Axx_(I), gf_Axy_(I), gf_Axz_(I),
             gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  template <typename F, typename = result_of_t<F(int, int)> >
  constexpr Z4C_INLINE mat3(const F &f)
      : elts{f(0, 0), f(0, 1), f(0, 2), f(1, 1), f(1, 2), f(2, 2)} {
    // #ifdef CCTK_DEBUG
    //     // Check symmetry
    //     const T f10 = f(1, 0);
    //     const T f20 = f(0, 2);
    //     const T f21 = f(1, 2);
    //     const auto is_symmetric{[](const T &fgood, const T &fother) {
    //       return norm1<T>()(fother - fgood) <=
    //              1.0e-12 * (1 + norm1<T>()(fgood) + norm1<T>()(fother));
    //     }};
    //     if (!(is_symmetric(f(0, 1), f10) && is_symmetric(f(0, 2), f20) &&
    //           is_symmetric(f(1, 2), f21))) {
    //       ostringstream buf;
    //       buf << "f(0,1)=" << f(0, 1) << "\n"
    //           << "f(1,0)=" << f10 << "\n"
    //           << "f(0,2)=" << f(0, 2) << "\n"
    //           << "f(2,0)=" << f20 << "\n"
    //           << "f(1,2)=" << f(1, 2) << "\n"
    //           << "f(2,1)=" << f21 << "\n";
    //       CCTK_VERROR("symmetric matrix is not symmetric:\n%s",
    //       buf.str().c_str());
    //     }
    //     assert(is_symmetric(f(0, 1), f10));
    //     assert(is_symmetric(f(0, 2), f20));
    //     assert(is_symmetric(f(1, 2), f21));
    // #endif
  }

  Z4C_INLINE void store(const GF3D2<T> &gf_Axx_, const GF3D2<T> &gf_Axy_,
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

  Z4C_INLINE const T &operator()(int i, int j) const { return elts[ind(i, j)]; }
  // Z4C_INLINE T &operator()(int i, int j) { return
  // elts[symind(i, j)]; }
  Z4C_INLINE T &operator()(int i, int j) { return elts[ind(i, j)]; }

  template <typename U = T>
  Z4C_INLINE
      mat3<remove_cv_t<remove_reference_t<result_of_t<U(vect<int, 3>)> > >,
           dnup1, dnup2>
      operator()(const vect<int, 3> &I) const {
    return {elts[0](I), elts[1](I), elts[2](I),
            elts[3](I), elts[4](I), elts[5](I)};
  }
  // TODO: Only if T is GF3D5<U>
  template <typename U = T>
  Z4C_INLINE
      mat3<remove_cv_t<
               remove_reference_t<result_of_t<U(GF3D5layout, vect<int, 3>)> > >,
           dnup1, dnup2>
      operator()(const GF3D5layout &layout, const vect<int, 3> &I) const {
    return {elts[0](layout, I), elts[1](layout, I), elts[2](layout, I),
            elts[3](layout, I), elts[4](layout, I), elts[5](layout, I)};
  }

  friend constexpr Z4C_INLINE mat3<T, dnup1, dnup2>
  operator+(const mat3<T, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend constexpr Z4C_INLINE mat3<T, dnup1, dnup2>
  operator-(const mat3<T, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend constexpr Z4C_INLINE mat3<T, dnup1, dnup2>
  operator+(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend constexpr Z4C_INLINE mat3<T, dnup1, dnup2>
  operator-(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend constexpr Z4C_INLINE mat3<T, dnup1, dnup2>
  operator*(const T &a, const mat3<T, dnup1, dnup2> &x) {
    return {a * x.elts};
  }

  friend constexpr Z4C_INLINE bool operator==(const mat3<T, dnup1, dnup2> &x,
                                              const mat3<T, dnup1, dnup2> &y) {
    return equal_to<vect<T, 6> >()(x.elts, y.elts);
  }
  friend constexpr Z4C_INLINE bool operator!=(const mat3<T, dnup1, dnup2> &x,
                                              const mat3<T, dnup1, dnup2> &y) {
    return !(x == y);
  }

  constexpr Z4C_INLINE T maxabs() const { return elts.maxabs(); }

  friend struct norm1<mat3>;

  constexpr Z4C_INLINE T det() const {
    const auto &A = *this;
    return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) -
           A(1, 0) * (A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) +
           A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
  }

  constexpr Z4C_INLINE mat3<T, !dnup1, !dnup2> inv(const T detA) const {
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

  constexpr Z4C_INLINE T trace(const mat3<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    return sum2([&](int x, int y) Z4C_INLINE { return gu(x, y) * A(x, y); });
  }

  constexpr Z4C_INLINE mat3 trace_free(
      const mat3<T, dnup1, dnup2> &g, const mat3<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    const T trA = A.trace(gu);
    return mat3([&](int a, int b)
                    Z4C_INLINE { return A(a, b) - trA / 3 * g(a, b); });
  }

  friend ostream &operator<<(ostream &os, const mat3<T, dnup1, dnup2> &A) {
    return os << "[[" << A(0, 0) << "," << A(0, 1) << "," << A(0, 2) << "],["
              << A(1, 0) << "," << A(1, 1) << "," << A(1, 2) << "],[" << A(2, 0)
              << "," << A(2, 1) << "," << A(2, 2) << "]]";
  }
};
template <typename T, dnup_t dnup1, dnup_t dnup2>
struct nan<mat3<T, dnup1, dnup2> > {
  constexpr Z4C_INLINE mat3<T, dnup1, dnup2> operator()() const {
    return mat3<T, dnup1, dnup2>();
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct norm1<mat3<T, dnup1, dnup2> > {
  typedef typename norm1<vect<T, 6> >::result_type result_type;
  constexpr Z4C_INLINE result_type
  operator()(const mat3<T, dnup1, dnup2> &x) const {
    return norm1<vect<T, 6> >()(x.elts);
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup3>
constexpr Z4C_INLINE mat3<T, dnup1, dnup2>
mul(const mat3<T, dnup1, dnup3> &A, const mat3<T, !dnup3, dnup2> &B) {
  // C[a,b] = A[a,c] B[c,b]
  return mat3<T, dnup1, dnup2>([&](int a, int b) Z4C_INLINE {
    return sum1([&](int x) Z4C_INLINE { return A(a, x) * B(x, b); });
  });
}
} // namespace Z4c

#endif // #ifndef TENSOR_HXX
