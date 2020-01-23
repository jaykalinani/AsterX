#ifndef TENSOR_HXX
#define TENSOR_HXX

#include <loop.hxx>

#include <cctk.h>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <ostream>
#include <sstream>
#include <type_traits>

namespace Z4c {
using namespace Loop;
using namespace std;

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T pow2(const T x) {
  return x * x;
}
template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T pow3(const T x) {
  return pow2(x) * x;
}
template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T pow4(const T x) {
  return pow2(pow2(x));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct nan {
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T operator()() const {
    return NAN;
  }
};
template <typename T, int D> struct nan<vect<T, D> > {
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vect<T, D> operator()() const {
    return vect<T, D>::pure(nan<T>()());
  }
};

template <typename T> struct norm1 {
  typedef T result_type;
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr result_type
  operator()(const T &x) const {
    return abs(x);
  }
};
template <typename T, int D> struct norm1<vect<T, D> > {
  typedef typename norm1<T>::result_type result_type;
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr result_type
  operator()(const vect<T, D> &xs) const {
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

  static /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr int ind(const int n) {
#ifdef CCTK_DEBUG
    assert(n >= 0 && n < 3);
#endif
    return n;
  }

public:
  explicit /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3()
      : elts{nan<vect<T, 3> >()()} {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3(const vect<T, 3> &elts)
      : elts(elts) {}
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3(vect<T, 3> &&elts)
      : elts(move(elts)) {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3(initializer_list<T> v)
      : elts(v) {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ vec3(const GF3D<const T, 0, 0, 0> &gf_vx_,
                                        const GF3D<const T, 0, 0, 0> &gf_vy_,
                                        const GF3D<const T, 0, 0, 0> &gf_vz_,
                                        const vect<int, 3> &I)
      : vec3{gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/
  vec3(const GF3D<T, 0, 0, 0> &gf_vx_, const GF3D<T, 0, 0, 0> &gf_vy_,
       const GF3D<T, 0, 0, 0> &gf_vz_, const vect<int, 3> &I)
      : vec3{gf_vx_(I), gf_vy_(I), gf_vz_(I)} {}

  template <typename F, typename = result_of_t<F(int)> >
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3(const F &f)
      : elts{f(0), f(1), f(2)} {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ void store(const GF3D<T, 0, 0, 0> &gf_vx_,
                                              const GF3D<T, 0, 0, 0> &gf_vy_,
                                              const GF3D<T, 0, 0, 0> &gf_vz_,
                                              const PointDesc &p) const {
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
    gf_vx_(p.I) = v(0);
    gf_vy_(p.I) = v(1);
    gf_vz_(p.I) = v(2);
  }

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ const T &operator()(int i) const {
    return elts[ind(i)];
  }
  T &operator()(int i) { return elts[ind(i)]; }

  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3<T, dnup>
  operator+(const vec3<T, dnup> &x) {
    return {+x.elts};
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3<T, dnup>
  operator-(const vec3<T, dnup> &x) {
    return {-x.elts};
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3<T, dnup>
  operator+(const vec3<T, dnup> &x, const vec3<T, dnup> &y) {
    return {x.elts + y.elts};
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3<T, dnup>
  operator-(const vec3<T, dnup> &x, const vec3<T, dnup> &y) {
    return {x.elts - y.elts};
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3<T, dnup>
  operator*(const T &a, const vec3<T, dnup> &x) {
    return {a * x.elts};
  }

  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr bool
  operator==(const vec3<T, dnup> &x, const vec3<T, dnup> &y) {
    return equal_to<vect<T, 3> >()(x.elts, y.elts);
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr bool
  operator!=(const vec3<T, dnup> &x, const vec3<T, dnup> &y) {
    return !(x == y);
  }

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T maxabs() const {
    return elts.maxabs();
  }

  friend struct norm1<vec3>;

  friend ostream &operator<<(ostream &os, const vec3<T, dnup> &v) {
    return os << "[" << v(0) << "," << v(1) << "," << v(2) << "]";
  }
};

template <typename T, dnup_t dnup> struct nan<vec3<T, dnup> > {
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr vec3<T, dnup> operator()() const {
    return vec3<T, dnup>();
  }
};

template <typename T, dnup_t dnup> struct norm1<vec3<T, dnup> > {
  typedef typename norm1<vect<T, 3> >::result_type result_type;
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr result_type
  operator()(const vec3<T, dnup> &x) const {
    return norm1<vect<T, 3> >()(x.elts);
  }
};

////////////////////////////////////////////////////////////////////////////////

// Symmetric 3-matrix
template <typename T, dnup_t dnup1, dnup_t dnup2> class mat3 {
  static_assert(dnup1 == dnup2, "");

  vect<T, 6> elts;

  static /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr int symind(const int i,
                                                               const int j) {
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
  static /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr int ind(const int i,
                                                            const int j) {
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
  explicit /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3()
      : elts{nan<vect<T, 6> >()()} {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3(const vect<T, 6> &elts)
      : elts(elts) {}
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3(vect<T, 6> &&elts)
      : elts(move(elts)) {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3(initializer_list<T> A)
      : elts(A) {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ mat3(const GF3D<const T, 0, 0, 0> &gf_Axx_,
                                        const GF3D<const T, 0, 0, 0> &gf_Axy_,
                                        const GF3D<const T, 0, 0, 0> &gf_Axz_,
                                        const GF3D<const T, 0, 0, 0> &gf_Ayy_,
                                        const GF3D<const T, 0, 0, 0> &gf_Ayz_,
                                        const GF3D<const T, 0, 0, 0> &gf_Azz_,
                                        const vect<int, 3> &I)
      : mat3{gf_Axx_(I), gf_Axy_(I), gf_Axz_(I),
             gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/
  mat3(const GF3D<T, 0, 0, 0> &gf_Axx_, const GF3D<T, 0, 0, 0> &gf_Axy_,
       const GF3D<T, 0, 0, 0> &gf_Axz_, const GF3D<T, 0, 0, 0> &gf_Ayy_,
       const GF3D<T, 0, 0, 0> &gf_Ayz_, const GF3D<T, 0, 0, 0> &gf_Azz_,
       const vect<int, 3> &I)
      : mat3{gf_Axx_(I), gf_Axy_(I), gf_Axz_(I),
             gf_Ayy_(I), gf_Ayz_(I), gf_Azz_(I)} {}

  template <typename F, typename = result_of_t<F(int, int)> >
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3(const F &f)
      : elts{f(0, 0), f(0, 1), f(0, 2), f(1, 1), f(1, 2), f(2, 2)} {
#ifdef CCTK_DEBUG
    // Check symmetry
    const T f10 = f(1, 0);
    const T f20 = f(0, 2);
    const T f21 = f(1, 2);
    const auto is_symmetric{[](const T &fgood, const T &fother) {
      return norm1<T>()(fother - fgood) <=
             1.0e-12 * (1 + norm1<T>()(fgood) + norm1<T>()(fother));
    }};
    if (!(is_symmetric(f(0, 1), f10) && is_symmetric(f(0, 2), f20) &&
          is_symmetric(f(1, 2), f21))) {
      ostringstream buf;
      buf << "f(0,1)=" << f(0, 1) << "\n"
          << "f(1,0)=" << f10 << "\n"
          << "f(0,2)=" << f(0, 2) << "\n"
          << "f(2,0)=" << f20 << "\n"
          << "f(1,2)=" << f(1, 2) << "\n"
          << "f(2,1)=" << f21 << "\n";
      CCTK_VERROR("symmetric matrix is not symmetric:\n%s", buf.str().c_str());
    }
    assert(is_symmetric(f(0, 1), f10));
    assert(is_symmetric(f(0, 2), f20));
    assert(is_symmetric(f(1, 2), f21));
#endif
  }

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ void
  store(const GF3D<T, 0, 0, 0> &gf_Axx_, const GF3D<T, 0, 0, 0> &gf_Axy_,
        const GF3D<T, 0, 0, 0> &gf_Axz_, const GF3D<T, 0, 0, 0> &gf_Ayy_,
        const GF3D<T, 0, 0, 0> &gf_Ayz_, const GF3D<T, 0, 0, 0> &gf_Azz_,
        const PointDesc &p) const {
    const auto &A = *this;
#ifdef CCTK_DEBUG
    assert(CCTK_isfinite(A(0, 0)));
    assert(CCTK_isfinite(A(0, 1)));
    assert(CCTK_isfinite(A(0, 2)));
    assert(CCTK_isfinite(A(1, 1)));
    assert(CCTK_isfinite(A(1, 2)));
    assert(CCTK_isfinite(A(2, 2)));
#endif
    gf_Axx_(p.I) = A(0, 0);
    gf_Axy_(p.I) = A(0, 1);
    gf_Axz_(p.I) = A(0, 2);
    gf_Ayy_(p.I) = A(1, 1);
    gf_Ayz_(p.I) = A(1, 2);
    gf_Azz_(p.I) = A(2, 2);
  }

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ const T &operator()(int i, int j) const {
    return elts[ind(i, j)];
  }
  // /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ T &operator()(int i, int j) { return
  // elts[symind(i, j)]; }
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ T &operator()(int i, int j) {
    return elts[ind(i, j)];
  }

  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<T, dnup1, dnup2>
  operator+(const mat3<T, dnup1, dnup2> &x) {
    return {+x.elts};
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<T, dnup1, dnup2>
  operator-(const mat3<T, dnup1, dnup2> &x) {
    return {-x.elts};
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<T, dnup1, dnup2>
  operator+(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return {x.elts + y.elts};
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<T, dnup1, dnup2>
  operator-(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return {x.elts - y.elts};
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<T, dnup1, dnup2>
  operator*(const T &a, const mat3<T, dnup1, dnup2> &x) {
    return {a * x.elts};
  }

  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr bool
  operator==(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return equal_to<vect<T, 6> >()(x.elts, y.elts);
  }
  friend /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr bool
  operator!=(const mat3<T, dnup1, dnup2> &x, const mat3<T, dnup1, dnup2> &y) {
    return !(x == y);
  }

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T maxabs() const {
    return elts.maxabs();
  }

  friend struct norm1<mat3>;

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T det() const {
    const auto &A = *this;
    return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) -
           A(1, 0) * (A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) +
           A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
  }

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<T, !dnup1, !dnup2>
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

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T
  trace(const mat3<T, !dnup1, !dnup2> &gu) const {
    const auto &A = *this;
    return sum2([&](int x, int y) { return gu(x, y) * A(x, y); });
  }

  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3 trace_free(
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
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<T, dnup1, dnup2>
  operator()() const {
    return mat3<T, dnup1, dnup2>();
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2>
struct norm1<mat3<T, dnup1, dnup2> > {
  typedef typename norm1<vect<T, 6> >::result_type result_type;
  /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr result_type
  operator()(const mat3<T, dnup1, dnup2> &x) const {
    return norm1<vect<T, 6> >()(x.elts);
  }
};

template <typename T, dnup_t dnup1, dnup_t dnup2, dnup_t dnup3>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr mat3<T, dnup1, dnup2>
mul(const mat3<T, dnup1, dnup3> &A, const mat3<T, !dnup3, dnup2> &B) {
  // C[a,b] = A[a,c] B[c,b]
  return mat3<T, dnup1, dnup2>([&](int a, int b) {
    return sum1([&](int x) { return A(a, x) * B(x, b); });
  });
}

template <typename F, typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T fold(const F &f, const T &x) {
  return x;
}
template <typename F, typename T, typename... Ts>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T fold(const F &f, const T &x0,
                                                  const T &x1,
                                                  const Ts &... xs) {
  return fold(f, fold(f, x0, x1), xs...);
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T add() {
  return T(0);
}
// template <typename T>
// /*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T add(const T &x) {
//   return x;
// }
template <typename T, typename... Ts>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr T add(const T &x, const Ts &... xs) {
  return x + add(xs...);
}

template <typename F>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr result_of_t<F(int)>
sum1(const F &f) {
  result_of_t<F(int)> s{0};
  for (int x = 0; x < 3; ++x)
    s += f(x);
  return s;
}

template <typename F>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr result_of_t<F(int, int)>
sum2(const F &f) {
  result_of_t<F(int, int)> s{0};
  for (int x = 0; x < 3; ++x)
    for (int y = 0; y < 3; ++y)
      s += f(x, y);
  return s;
}

template <typename F>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ constexpr result_of_t<F(int, int, int)>
sum3(const F &f) {
  result_of_t<F(int, int, int)> s{0};
  for (int x = 0; x < 3; ++x)
    for (int y = 0; y < 3; ++y)
      for (int z = 0; z < 3; ++z)
        s += f(x, y, z);
  return s;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ T deriv(const GF3D<const T, 0, 0, 0> &gf_,
                                         const vect<int, dim> &I,
                                         const vec3<T, UP> &dx, const int dir) {
  constexpr vect<vect<int, dim>, dim> DI{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  return (gf_(I + DI[dir]) - gf_(I - DI[dir])) / (2 * dx(dir));
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ T
deriv_upwind(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
             const bool sign, const vec3<T, UP> &dx, const int dir) {
  constexpr vect<vect<int, dim>, dim> DI{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  // arXiv:1111.2177 [gr-qc], (71)
  if (sign)
    // +     [ 0   -1   +1    0    0]
    // + 1/2 [+1   -2   +1    0    0]
    //       [+1/2 -2   +3/2  0    0]
    return (+1 / T(2) * gf_(I - 2 * DI[dir]) //
            - 2 * gf_(I - DI[dir])           //
            + 3 / T(2) * gf_(I)) /
           dx(dir);
  else
    // +     [ 0    0   -1   +1    0  ]
    // - 1/2 [ 0    0   +1   -2   +1  ]
    //       [ 0    0   -3/2 +2   -1/2]
    return (-3 / T(2) * gf_(I)     //
            + 2 * gf_(I + DI[dir]) //
            - 1 / T(2) * gf_(I + 2 * DI[dir])) /
           dx(dir);
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ T
deriv2(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
       const vec3<T, UP> &dx, const int dir1, const int dir2) {
  constexpr vect<vect<int, dim>, dim> DI{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  if (dir1 == dir2)
    return ((gf_(I - DI[dir1]) + gf_(I + DI[dir1])) - 2 * gf_(I)) /
           pow2(dx(dir1));
  else
    // return (deriv(gf_, I + DI[dir2], dx, dir1) -
    //         deriv(gf_, I - DI[dir2], dx, dir1)) /
    //        (2 * dx(dir2));
    return ((gf_(I + DI[dir1] + DI[dir2])       //
             - gf_(I - DI[dir1] + DI[dir2]))    //
            - (gf_(I + DI[dir1] - DI[dir2])     //
               - gf_(I - DI[dir1] - DI[dir2]))) //
           / (4 * dx(dir1) * dx(dir2));
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ T
deriv4(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
       const vec3<T, UP> &dx, const int dir) {
  constexpr vect<vect<int, dim>, dim> DI{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  return ((gf_(I - 2 * DI[dir]) + gf_(I + 2 * DI[dir])) //
          - 4 * (gf_(I - DI[dir]) + gf_(I + DI[dir]))   //
          + 6 * gf_(I)) /
         pow4(dx(dir));
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ vec3<T, DN>
deriv(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
      const vec3<T, UP> &dx) {
  return {
      deriv(gf_, I, dx, 0),
      deriv(gf_, I, dx, 1),
      deriv(gf_, I, dx, 2),
  };
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ vec3<T, DN>
deriv_upwind(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
             const vec3<T, UP> &dir, const vec3<T, UP> &dx) {
  return {
      deriv_upwind(gf_, I, signbit(dir(0)), dx, 0),
      deriv_upwind(gf_, I, signbit(dir(1)), dx, 1),
      deriv_upwind(gf_, I, signbit(dir(2)), dx, 2),
  };
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ mat3<T, DN, DN>
deriv2(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
       const vec3<T, UP> &dx) {
  return {
      deriv2(gf_, I, dx, 0, 0), deriv2(gf_, I, dx, 0, 1),
      deriv2(gf_, I, dx, 0, 2), deriv2(gf_, I, dx, 1, 1),
      deriv2(gf_, I, dx, 1, 2), deriv2(gf_, I, dx, 2, 2),
  };
}

template <typename T>
/*CCTK_ATTRIBUTE_ALWAYS_INLINE*/ T diss(const GF3D<const T, 0, 0, 0> &gf_,
                                        const vect<int, dim> &I,
                                        const vec3<T, UP> &dx) {
  // arXiv:gr-qc/0610128, (63), with r=2
  return -1 / T(16) *
         (pow3(dx(0)) * deriv4(gf_, I, dx, 0)   //
          + pow3(dx(1)) * deriv4(gf_, I, dx, 1) //
          + pow3(dx(2)) * deriv4(gf_, I, dx, 2));
}
} // namespace Z4c

#endif // #ifndef TENSOR_HXX
