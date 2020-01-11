#ifndef TENSOR_HXX
#define TENSOR_HXX

#include <loop.hxx>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <ostream>

namespace Z4c {
using namespace Loop;
using namespace std;

template <typename T> constexpr T pow2(const T x) { return x * x; }

////////////////////////////////////////////////////////////////////////////////

template <typename T> struct nan {
  constexpr T operator()() const { return 0.0 / 0.0; }
};
template <typename T, int D> struct nan<vect<T, D> > {
  constexpr vect<T, D> operator()() const {
    return vect<T, D>::pure(nan<T>()());
  }
};

////////////////////////////////////////////////////////////////////////////////

// 3-vector
template <typename T> class vec3 {
  vect<T, 3> elts;

  static constexpr int ind(const int n) {
    assert(n >= 0 && n < 3);
    return n;
  }

public:
  explicit constexpr vec3() : elts{nan<vect<T, 3> >()()} {}

  constexpr vec3(initializer_list<T> v) : elts(v) {}

  vec3(const GF3D<const T, 0, 0, 0> &gf_vx_,
       const GF3D<const T, 0, 0, 0> &gf_vy_,
       const GF3D<const T, 0, 0, 0> &gf_vz_, const PointDesc &p)
      : vec3{gf_vx_(p.I), gf_vy_(p.I), gf_vz_(p.I)} {}

  vec3(const GF3D<T, 0, 0, 0> &gf_vx_, const GF3D<T, 0, 0, 0> &gf_vy_,
       const GF3D<T, 0, 0, 0> &gf_vz_, const PointDesc &p)
      : vec3{gf_vx_(p.I), gf_vy_(p.I), gf_vz_(p.I)} {}

  void store(const GF3D<T, 0, 0, 0> &gf_vx_, const GF3D<T, 0, 0, 0> &gf_vy_,
             const GF3D<T, 0, 0, 0> &gf_vz_, const PointDesc &p) const {
    const auto &v = *this;
#ifdef CCTK_DEBUG
    assert(!CCTK_isnan(v(0)));
    assert(!CCTK_isnan(v(1)));
    assert(!CCTK_isnan(v(2)));
#endif
    gf_vx_(p.I) = v(0);
    gf_vy_(p.I) = v(1);
    gf_vz_(p.I) = v(2);
  }

  const T &operator()(int i) const { return elts[ind(i)]; }
  T &operator()(int i) { return elts[ind(i)]; }

  friend ostream &operator<<(ostream &os, const vec3<T> &v) {
    return os << "[" << v(0) << "," << v(1) << "," << v(2) << "]";
  }
};
template <typename T> struct nan<vec3<T> > {
  constexpr vec3<T> operator()() const { return vec3<T>(); }
};

////////////////////////////////////////////////////////////////////////////////

// Symmetric 3-matrix
template <typename T> class mat3 {
  vect<T, 6> elts;

  static constexpr int symind(const int i, const int j) {
    assert(i >= 0 && i <= j && j < 3);
    const int n = 3 * i - i * (i + 1) / 2 + j;
    assert(n >= 0 && n < 6);
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

  static_assert(ind(1, 0) == ind(0, 1), "");
  static_assert(ind(2, 0) == ind(0, 2), "");
  static_assert(ind(2, 1) == ind(1, 2), "");

public:
  explicit constexpr mat3() : elts{nan<vect<T, 6> >()()} {}

  constexpr mat3(initializer_list<T> A) : elts(A) {}

  mat3(const GF3D<const T, 0, 0, 0> &gf_Axx_,
       const GF3D<const T, 0, 0, 0> &gf_Axy_,
       const GF3D<const T, 0, 0, 0> &gf_Axz_,
       const GF3D<const T, 0, 0, 0> &gf_Ayy_,
       const GF3D<const T, 0, 0, 0> &gf_Ayz_,
       const GF3D<const T, 0, 0, 0> &gf_Azz_, const PointDesc &p)
      : mat3{gf_Axx_(p.I), gf_Axy_(p.I), gf_Axz_(p.I),
             gf_Ayy_(p.I), gf_Ayz_(p.I), gf_Azz_(p.I)} {}

  mat3(const GF3D<T, 0, 0, 0> &gf_Axx_, const GF3D<T, 0, 0, 0> &gf_Axy_,
       const GF3D<T, 0, 0, 0> &gf_Axz_, const GF3D<T, 0, 0, 0> &gf_Ayy_,
       const GF3D<T, 0, 0, 0> &gf_Ayz_, const GF3D<T, 0, 0, 0> &gf_Azz_,
       const PointDesc &p)
      : mat3{gf_Axx_(p.I), gf_Axy_(p.I), gf_Axz_(p.I),
             gf_Ayy_(p.I), gf_Ayz_(p.I), gf_Azz_(p.I)} {}

  void store(const GF3D<T, 0, 0, 0> &gf_Axx_, const GF3D<T, 0, 0, 0> &gf_Axy_,
             const GF3D<T, 0, 0, 0> &gf_Axz_, const GF3D<T, 0, 0, 0> &gf_Ayy_,
             const GF3D<T, 0, 0, 0> &gf_Ayz_, const GF3D<T, 0, 0, 0> &gf_Azz_,
             const PointDesc &p) const {
    const auto &A = *this;
#ifdef CCTK_DEBUG
    assert(!CCTK_isnan(A(0, 0)));
    assert(!CCTK_isnan(A(0, 1)));
    assert(!CCTK_isnan(A(0, 2)));
    assert(!CCTK_isnan(A(1, 1)));
    assert(!CCTK_isnan(A(1, 2)));
    assert(!CCTK_isnan(A(2, 2)));
#endif
    gf_Axx_(p.I) = A(0, 0);
    gf_Axy_(p.I) = A(0, 1);
    gf_Axz_(p.I) = A(0, 2);
    gf_Ayy_(p.I) = A(1, 1);
    gf_Ayz_(p.I) = A(1, 2);
    gf_Azz_(p.I) = A(2, 2);
  }

  const T &operator()(int i, int j) const { return elts[ind(i, j)]; }
  // T &operator()(int i, int j) { return elts[symind(i, j)]; }
  T &operator()(int i, int j) { return elts[ind(i, j)]; }

  constexpr T det() const {
    const auto &A = *this;
    return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) -
           A(1, 0) * (A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1)) +
           A(2, 0) * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
  }

  constexpr mat3 inv(const T detA) const {
    const auto &A = *this;
    const T detA1 = 1 / detA;
    return mat3{detA1 * (A(1, 1) * A(2, 2) - A(2, 1) * A(2, 1)),
                detA1 * (A(1, 2) * A(2, 0) - A(2, 2) * A(0, 1)),
                detA1 * (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1)),
                detA1 * (A(2, 2) * A(0, 0) - A(0, 2) * A(2, 0)),
                detA1 * (A(2, 0) * A(0, 1) - A(0, 0) * A(2, 1)),
                detA1 * (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0))};
  }

  constexpr T trace() const {
    const auto &A = *this;
    return A(0, 0) + A(1, 1) + A(2, 2);
  }

  friend ostream &operator<<(ostream &os, const mat3<T> &A) {
    return os << "[[" << A(0, 0) << "," << A(0, 1) << "," << A(0, 2) << "],["
              << A(1, 0) << "," << A(1, 1) << "," << A(1, 2) << "],[" << A(2, 0)
              << "," << A(2, 1) << "," << A(2, 2) << "]]";
  }
};
template <typename T> struct nan<mat3<T> > {
  constexpr mat3<T> operator()() const { return mat3<T>(); }
};

template <typename T>
constexpr mat3<T> mul(const mat3<T> &A, const mat3<T> &B) {
  mat3<T> C;
  // C[a,b] = A[a,c] B[c,b]
  for (int a = 0; a < 3; ++a)
    for (int b = a; b < 3; ++b) {
      T s = 0;
      for (int c = 0; c < 3; ++c)
        s += A(a, c) * B(c, b);
      C(a, b) = s;
    }
  return C;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
T deriv(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
        const vec3<T> &dx, const int dir) {
  constexpr vect<vect<int, dim>, dim> DI{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  return (gf_(I + DI[dir]) - gf_(I - DI[dir])) / (2 * dx(dir));
}

template <typename T>
vec3<T> deriv(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
              const vec3<T> &dx) {
  return {
      deriv(gf_, I, dx, 0),
      deriv(gf_, I, dx, 1),
      deriv(gf_, I, dx, 2),
  };
}

template <typename T>
T deriv2(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
         const vec3<T> &dx, const int dir1, const int dir2) {
  constexpr vect<vect<int, dim>, dim> DI{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  if (dir1 == dir2)
    return (gf_(I + DI[dir1]) - 2 * gf_(I) + gf_(I - DI[dir1])) /
           (dx(dir1) * dx(dir2));
  else
    return (deriv(gf_, I + DI[dir2], dx, dir1) -
            deriv(gf_, I - DI[dir2], dx, dir1)) /
           (2 * dx(dir2));
}

template <typename T>
mat3<T> deriv2(const GF3D<const T, 0, 0, 0> &gf_, const vect<int, dim> &I,
               const vec3<T> &dx) {
  return {
      deriv2(gf_, I, dx, 0, 0), deriv2(gf_, I, dx, 0, 1),
      deriv2(gf_, I, dx, 0, 2), deriv2(gf_, I, dx, 1, 1),
      deriv2(gf_, I, dx, 1, 2), deriv2(gf_, I, dx, 2, 2),
  };
}

} // namespace Z4c

#endif // #ifndef TENSOR_HXX
