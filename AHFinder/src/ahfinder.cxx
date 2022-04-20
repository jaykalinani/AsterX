#include "discretization.hxx"
#include "physics.hxx"

#include <dual.hxx>
#include <mat.hxx>
#include <sum.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <ssht/ssht.h>

#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

namespace AHFinder {
using namespace Arith;

template <typename T> using vec3 = vec<T, 3, UP>;
template <typename T> using mat3 = mat<T, 3, DN, DN>;
template <typename T> using smat3 = smat<T, 3, DN, DN>;

template <typename F> auto sum3(F &&f) { return sum<3>(std::forward<F>(f)); }

template <typename T> struct coords_t {
  geom_t geom;
  vec3<scalar_aij_t<T> > x;

  coords_t() = delete;
  coords_t(const geom_t &geom)
      : geom(geom), x{
                        scalar_aij_t<T>(geom),
                        scalar_aij_t<T>(geom),
                        scalar_aij_t<T>(geom),
                    } {}
};

template <typename T> coords_t<T> coords_from_shape(const scalar_aij_t<T> &h) {
  DECLARE_CCTK_PARAMETERS;
  using std::cos, std::sin;
  const geom_t &geom = h.geom;
  coords_t<T> coords(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      const T r = h()(i, j);
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      coords.x(0)()(i, j) = x0 + r * sin(theta) * cos(phi);
      coords.x(1)()(i, j) = y0 + r * sin(theta) * sin(phi);
      coords.x(2)()(i, j) = z0 + r * cos(theta);
    }
  }
  return coords;
}

template <typename T> struct metric_t {
  geom_t geom;
  smat3<scalar_aij_t<T> > g, K;
  smat3<vec3<scalar_aij_t<T> > > dg;

  metric_t() = delete;
  metric_t(const geom_t &geom)
      : geom(geom),
        g{
            scalar_aij_t<T>(geom), scalar_aij_t<T>(geom), scalar_aij_t<T>(geom),
            scalar_aij_t<T>(geom), scalar_aij_t<T>(geom), scalar_aij_t<T>(geom),
        },
        K{
            scalar_aij_t<T>(geom), scalar_aij_t<T>(geom), scalar_aij_t<T>(geom),
            scalar_aij_t<T>(geom), scalar_aij_t<T>(geom), scalar_aij_t<T>(geom),
        },
        dg{
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
            {
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
                scalar_aij_t<T>(geom),
            },
        } {}
};

template <typename T>
metric_t<T> brill_lindquist_metric(const cGH *const cctkGH,
                                   const coords_t<T> &coords) {
  const geom_t &geom = coords.geom;
  metric_t<T> metric(geom);

  const T M = 1.0;
  const vec3<T> x0{0.0, 0.0, 0.0};
  const mat3<T> fx{
      1.0, 0.0, 0.0, //
      0.0, 1.0, 0.0, //
      0.0, 0.0, 1.0, //
  };

  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      if (0)
        std::cout << "i:" << i << " j:" << j << "\n";
      using U = vect<T, 3>;
      using DT = dual<T, U>; // spatial derivatives (x, y, z)
      const vec3<DT> x{
          {coords.x(0)()(i, j), {1, 0, 0}},
          {coords.x(1)()(i, j), {0, 1, 0}},
          {coords.x(2)()(i, j), {0, 0, 1}},
      };
      if (0)
        std::cout << "  x:" << x << "\n";
      const vec3<DT> X = matmul(fx, x - x0);
      if (0)
        std::cout << "  X:" << X << "\n";
      using std::pow, std::sqrt;
      const DT R = sqrt(pow(X(0), 2) + pow(X(1), 2) + pow(X(2), 2));
      const DT Psi = pow(1 + M / (2 * R), 4);
      const smat3<T> eta([&](int a, int b) { return pow(X(a).eps[b], 2); });
      if (0)
        std::cout << "  eta:" << eta << "\n";
      const smat3<DT> gg([&](int a, int b) { return Psi * eta(a, b); });
      const smat3<T> g([&](int a, int b) { return gg(a, b).val; });
      const smat3<vec3<T> > dg([&](int a, int b) {
        return vec3<T>([&](int c) { return gg(a, b).eps[c]; });
      });
      if (0)
        std::cout << "  g:" << g << "\n";
      if (0)
        std::cout << "  dg:" << dg << "\n";
      const smat3<T> K([&](int a, int b) { return 0; });
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          metric.g(a, b)()(i, j) = g(a, b);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          for (int c = 0; c < 3; ++c)
            metric.dg(a, b)(c)()(i, j) = dg(a, b)(c);
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b)
          metric.K(a, b)()(i, j) = K(a, b);
    }
  }

  return metric;
}

template <typename T>
metric_t<T> interpolate_metric(const cGH *const cctkGH,
                               const coords_t<T> &coords) {
  const int gxx_ind = CCTK_VarIndex("ADMBase::gxx");
  const int gxy_ind = CCTK_VarIndex("ADMBase::gxy");
  const int gxz_ind = CCTK_VarIndex("ADMBase::gxz");
  const int gyy_ind = CCTK_VarIndex("ADMBase::gyy");
  const int gyz_ind = CCTK_VarIndex("ADMBase::gyz");
  const int gzz_ind = CCTK_VarIndex("ADMBase::gzz");
  const int kxx_ind = CCTK_VarIndex("ADMBase::kxx");
  const int kxy_ind = CCTK_VarIndex("ADMBase::kxy");
  const int kxz_ind = CCTK_VarIndex("ADMBase::kxz");
  const int kyy_ind = CCTK_VarIndex("ADMBase::kyy");
  const int kyz_ind = CCTK_VarIndex("ADMBase::kyz");
  const int kzz_ind = CCTK_VarIndex("ADMBase::kzz");

  constexpr int nvars = 6 * (1 + 3 + 1);
  const array<CCTK_INT, nvars> varinds{
      gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
      gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
      gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
      gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
      kxx_ind, kxy_ind, kxz_ind, kyy_ind, kyz_ind, kzz_ind, //
  };
  const array<CCTK_INT, nvars> operations{
      0, 0, 0, 0, 0, 0, //
      1, 1, 1, 1, 1, 1, //
      2, 2, 2, 2, 2, 2, //
      3, 3, 3, 3, 3, 3, //
      0, 0, 0, 0, 0, 0, //
  };

  const geom_t &geom = coords.geom;
  metric_t<T> metric(geom);
  array<T *, nvars> ptrs{
      metric.g(0, 0)().data(),     metric.g(0, 1)().data(),
      metric.g(0, 2)().data(),     metric.g(1, 1)().data(),
      metric.g(1, 2)().data(),     metric.g(2, 2)().data(),
      metric.dg(0, 0)(0)().data(), metric.dg(0, 1)(0)().data(),
      metric.dg(0, 2)(0)().data(), metric.dg(1, 1)(0)().data(),
      metric.dg(1, 2)(0)().data(), metric.dg(2, 2)(0)().data(),
      metric.dg(0, 0)(1)().data(), metric.dg(0, 1)(1)().data(),
      metric.dg(0, 2)(1)().data(), metric.dg(1, 1)(1)().data(),
      metric.dg(1, 2)(1)().data(), metric.dg(2, 2)(1)().data(),
      metric.dg(0, 0)(2)().data(), metric.dg(0, 1)(2)().data(),
      metric.dg(0, 2)(2)().data(), metric.dg(1, 1)(2)().data(),
      metric.dg(1, 2)(2)().data(), metric.dg(2, 2)(2)().data(),
      metric.K(0, 0)().data(),     metric.K(0, 1)().data(),
      metric.K(0, 2)().data(),     metric.K(1, 1)().data(),
      metric.K(1, 2)().data(),     metric.K(2, 2)().data()};

  Interpolate(cctkGH, geom.npoints, coords.x(0)().data(), coords.x(1)().data(),
              coords.x(2)().data(), nvars, varinds.data(), operations.data(),
              ptrs.data());

  metric.g(1, 0)() = metric.g(0, 1)();
  metric.g(2, 0)() = metric.g(0, 2)();
  metric.g(2, 1)() = metric.g(1, 2)();
  metric.dg(1, 0)() = metric.dg(0, 1)();
  metric.dg(2, 0)() = metric.dg(0, 2)();
  metric.dg(2, 1)() = metric.dg(1, 2)();
  metric.K(1, 0)() = metric.K(0, 1)();
  metric.K(2, 0)() = metric.K(0, 2)();
  metric.K(2, 1)() = metric.K(1, 2)();

  return metric;
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
scalar_alm_t<std::complex<T> >
expansion(const metric_t<T> &metric,
          const scalar_alm_t<std::complex<T> > &hlm) {
  DECLARE_CCTK_PARAMETERS;

  using std::cos, std::sin;

  const geom_t &geom = hlm.geom;

  // Cartesian coordinates:
  //     x = r \cos\theta \cos\phi
  //     y = r \cos\theta \sin\phi
  //     z = r \sin\theta

  // Surface coordinates:
  //     F(r,\theta,\phi) = r - h(\theta,\phi)
  //     F   = r - h(\theta,\phi)
  //     \mu = \theta
  //     \nu = \phi

  // Derivatives Cartesian/spherical:
  //     dx/dr      =     \cos\theta \cos\phi
  //     dy/dr      =     \cos\theta \sin\phi
  //     dz/dr      =     \sin\theta
  //     dx/d\theta = - r \sin\theta \cos\phi
  //     dy/d\theta = - r \sin\theta \sin\phi
  //     dz/d\theta =   r \cos\theta
  //     dx/d\phi   = - r \cos\theta \sin\phi
  //     dy/d\phi   =   r \cos\theta \cos\phi
  //     dz/d\phi   =   0

  // Derivatives surface/spherical:
  //     dF  /dr      = 1
  //     d\mu/dr      = 0
  //     d\nu/dr      = 0
  //     df  /d\theta = - dh/d\theta(\theta,\phi)
  //     d\mu/d\theta = 1
  //     d\nu/d\theta = 0
  //     df  /d\phi   = - dh/d\phi(\theta,\phi)
  //     d\mu/d\phi   = 0
  //     d\nu/d\phi   = 1

  // Evaluate h and its derivatives
  const scalar_aij_t<T> hij = evaluate(hlm);
  const vector_alm_t<std::complex<T> > dhlm = gradient(hlm);
  const vector_aij_t<T> dhij = evaluate(dhlm);
  if (0)
    std::cout << "dhij:" << dhij << "\n";

  // s^a
  const vec3<scalar_aij_t<T> > sij{
      scalar_aij_t<T>(geom),
      scalar_aij_t<T>(geom),
      scalar_aij_t<T>(geom),
  };
  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      if (0)
        std::cout << "i:" << i << " j:" << j << "\n";

      // Coordinates
      const T r = hij()(i, j);
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      if (0)
        std::cout << "  r:" << r << " theta:" << theta << " phi:" << phi
                  << "\n";
      const vec3<T> x{
          x0 + r * sin(theta) * cos(phi),
          y0 + r * sin(theta) * sin(phi),
          z0 + r * cos(theta),
      };
      if (0)
        std::cout << "  x:" << x << "\n";
      const mat3<T> dxdr{
          sin(theta) * cos(phi),
          r * cos(theta) * cos(phi),
          -r * sin(phi), // dx/dphi / sin(theta)

          sin(theta) * sin(phi),
          r * cos(theta) * sin(phi),
          r * cos(phi), // dz/dphi / sin(theta)

          cos(theta),
          -r * sin(theta),
          0, // dz/dphi / sin(theta)
      };
      if (0)
        std::cout << "  dxdr:" << dxdr << "\n";
      const mat3<T> drdx = inv(dxdr);
      if (0)
        std::cout << "  drdx:" << drdx << "\n";

      // Level set function
      // F(r, theta, phi) = r - h(theta, phi)
      // const T F = 0
      const mat3<T> dFdr{
          1,              // dF/dr
          -dhij(0)(i, j), // dF/dtheta
          -dhij(1)(i, j), // dF/dphi / sin(theta)

          0,
          1,
          0,

          0,
          0,
          1,
      };
      if (0)
        std::cout << "  dFdr:" << dFdr << "\n";
      const mat3<T> dFdx([&](int a, int b) {
        return sum3([&](int x) { return dFdr(a, x) * drdx(x, b); });
      });

      // Metric
      const smat3<T> g([&](int a, int b) { return metric.g(a, b)()(i, j); });
      const smat3<T> gu = inv(g);

      // Surface normal
      const vec3<T> grad_F_F{1, 0, 0};

      const vec3<T> grad_F([&](int a) {
        return sum3([&](int x) { return grad_F_F(x) * dFdx(x, a); });
      });
      const vec3<T> grad_F_u([&](int a) {
        return sum3([&](int x) { return gu(a, x) * grad_F(x); });
      });

      const T len2_grad_F =
          sum3([&](int x) { return grad_F_u(x) * grad_F(x); });
      const vec3<T> s = grad_F_u / sqrt(len2_grad_F);
      for (int a = 0; a < 3; ++a)
        (T &)sij(a)()(i, j) = s(a);
    }
  }

  const vec3<scalar_alm_t<std::complex<T> > > slm(
      [&](int a) { return expand(sij(a)); });
  const vec3<vector_alm_t<std::complex<T> > > dslm(
      [&](int a) { return gradient(slm(a)); });
  const vec3<vector_aij_t<T> > dsij([&](int a) { return evaluate(dslm(a)); });
  if (0)
    std::cout << "dsij:" << dsij << "\n";

  scalar_aij_t<T> Thetaij(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      if (0)
        std::cout << "i:" << i << " j:" << j << "\n";

      // Coordinates
      const T r = hij()(i, j);
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      if (0)
        std::cout << "  r:" << r << " theta:" << theta << " phi:" << phi
                  << "\n";
      const vec3<T> x{
          x0 + r * sin(theta) * cos(phi),
          y0 + r * sin(theta) * sin(phi),
          z0 + r * cos(theta),
      };
      if (0)
        std::cout << "  x:" << x << "\n";
      const mat3<T> dxdr{
          sin(theta) * cos(phi),
          r * cos(theta) * cos(phi),
          -r * sin(phi), // dx/dphi / sin(theta)

          sin(theta) * sin(phi),
          r * cos(theta) * sin(phi),
          r * cos(phi), // dz/dphi / sin(theta)

          cos(theta),
          -r * sin(theta),
          0, // dz/dphi / sin(theta)
      };
      if (0)
        std::cout << "  dxdr:" << dxdr << "\n";
      const mat3<T> drdx = inv(dxdr);
      if (0)
        std::cout << "  drdx:" << drdx << "\n";

      // Level set function
      // F(r, theta, phi) = r - h(theta, phi)
      // const T F = 0
      const mat3<T> dFdr{
          1,              // dF/dr
          -dhij(0)(i, j), // dF/dtheta
          -dhij(1)(i, j), // dF/dphi / sin(theta)

          0,
          1,
          0,

          0,
          0,
          1,
      };
      if (0)
        std::cout << "  dFdr:" << dFdr << "\n";
      const mat3<T> dFdx([&](int a, int b) {
        return sum3([&](int x) { return dFdr(a, x) * drdx(x, b); });
      });

      // Metric
      const smat3<T> g([&](int a, int b) { return metric.g(a, b)()(i, j); });
      const smat3<T> K([&](int a, int b) { return metric.K(a, b)()(i, j); });
      const smat3<T> gu = inv(g);
      const smat3<vec3<T> > dg([&](int a, int b) {
        return vec3<T>([&](int c) { return metric.dg(a, b)(c)()(i, j); });
      });
      const vec3<smat3<T> > Gamma([&](int a) {
        return smat3<T>([&](int b, int c) {
          return sum3([&](int x) {
            return gu(a, x) * (dg(x, c)(b) + dg(b, x)(c) - dg(b, c)(x)) / 2;
          });
        });
      });

      // Surface normal
      const vec3<T> s([&](int a) { return sij(a)()(i, j); });
      if (0)
        std::cout << "  s:" << s << "\n";
      const mat3<T> dsdF{
          0, dsij(0)(0)(i, j), dsij(0)(1)(i, j),

          0, dsij(1)(0)(i, j), dsij(1)(1)(i, j),

          0, dsij(2)(0)(i, j), dsij(2)(1)(i, j),
      };
      if (0)
        std::cout << "  dsdF:" << dsdF << "\n";
      const mat3<T> dsdx([&](int a, int b) {
        return sum3([&](int x) { return dsdF(a, x) * dFdx(x, b); });
      });
      if (0)
        std::cout << "  dsdx:" << dsdx << "\n";
      const mat3<T> grad_s([&](int a, int b) {
        return dsdx(a, b) + sum3([&](int x) { return Gamma(a)(b, x) * s(x); });
      });

      // Expansion Theta_(l)
      // Change sign in front of `K` for Theta_(n)
      const T Theta = sum3([&](int x, int y) {
        return (gu(x, y) - s(x) * s(y)) * (grad_s(x, y) - K(x, y));
      });
      if (0)
        std::cout << "  Theta:" << Theta << "\n";
      Thetaij()(i, j) = Theta;
    }
  }

  const scalar_alm_t<std::complex<T> > Thetalm = expand(Thetaij);

  return Thetalm;

#if 0

  const T area = sqrt(4 * T(M_PI)) * real(Alm(0, 0));

  const T cx = real(cxlm(0, 0)) / real(Alm(0, 0));
  const T cy = real(cylm(0, 0)) / real(Alm(0, 0));
  const T cz = real(czlm(0, 0)) / real(Alm(0, 0));

  const alm_t<std::complex<T> > Thetalm = expand(Thetaij, 0);
  // [arXiv:gr-qc/0702038], (28)
  aij_t<T> Sij(geom);
  for (int i = 0; i < geom.ntheta; ++i)
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j)
      Sij(i, j) = lambdaij(i, j) * Thetaij(i, j);
  const alm_t<std::complex<T> > Slm = expand(Sij, 0);

  alm_t<std::complex<T> > hlm_new(geom, 0);
  for (int l = 0; l <= geom.lmax; ++l) {
#pragma omp simd
    for (int m = -l; m <= l; ++m) {
      hlm_new(l, m) = hlm(l, m) - 1 / T(l * (l + 1) + 2) * Slm(l, m);
    }
  }

  return {.hlm = hlm,
          .area = area,
          .cx = cx,
          .cy = cy,
          .cz = cz,
          .Thetalm = Thetalm,
          .hlm_new = hlm_new};

#endif
}

////////////////////////////////////////////////////////////////////////////////

#if 0
  
template <typename T>
expansion_t<T> update(const cGH *const cctkGH,
                      const alm_t<std::complex<T> > &hlm) {
  const auto hij = evaluate(hlm);
  const auto coords = coords_from_shape(hij);
  // CCTK_VWARN(CCTK_WARN_ALERT, "Using Brill-Lindquist metric");
  // const auto metric = brill_lindquist_metric(cctkGH, coords);
  const auto metric = interpolate_metric(cctkGH, coords);
  const auto res = expansion(metric, hlm);
  return res;
}

template <typename T>
expansion_t<T> solve(const cGH *const cctkGH,
                     const alm_t<std::complex<T> > &hlm_ini) {
  DECLARE_CCTK_PARAMETERS;
  int iter = 0;
  unique_ptr<const alm_t<std::complex<T> > > hlm_ptr =
      make_unique<alm_t<std::complex<T> > >(filter(hlm_ini, lmax_filter));
  for (;;) {
    ++iter;
    const auto &hlm = *hlm_ptr;
    const geom_t &geom = hlm.geom;
    const auto hij = evaluate(hlm);

    const auto res = update(cctkGH, hlm);

    const auto &Thetalm = res.Thetalm;
    const auto hlm_new = filter(res.hlm_new, lmax_filter);
    const auto hij_new = evaluate(hlm_new);

    T dh_maxabs{0};
    for (int i = 0; i < geom.ntheta; ++i)
#pragma omp simd
      for (int j = 0; j < geom.nphi; ++j)
        dh_maxabs = fmax(dh_maxabs, fabs(hij_new(i, j) - hij(i, j)));

    alm_t<std::complex<T> > dhlm(geom, 0);
    for (int l = 0; l <= geom.lmax; ++l)
#pragma omp simd
      for (int m = -l; m <= l; ++m)
        dhlm(l, m) = hlm_new(l, m) - hlm(l, m);
    auto dhij = evaluate(dhlm);
    aij_t<T> dh2ij(geom);
    for (int i = 0; i < geom.ntheta; ++i)
#pragma omp simd
      for (int j = 0; j < geom.nphi; ++j)
        dh2ij(i, j) = pow(dhij(i, j), 2);
    const auto dh2lm = expand(dh2ij, 0);
    const T dh_norm2 = sqrt(sqrt(4 * M_PI) * real(dh2lm(0, 0)));

    const auto Thetaij = evaluate(Thetalm);
    T Theta_maxabs{0};
    for (int i = 0; i < geom.ntheta; ++i)
#pragma omp simd
      for (int j = 0; j < geom.nphi; ++j)
        Theta_maxabs = fmax(Theta_maxabs, fabs(Thetaij(i, j)));

    aij_t<T> Theta2ij(geom);
    for (int i = 0; i < geom.ntheta; ++i)
#pragma omp simd
      for (int j = 0; j < geom.nphi; ++j)
        Theta2ij(i, j) = pow(Thetaij(i, j), 2);
    const auto Theta2lm = expand(Theta2ij, 0);
    const T Theta_norm2 = sqrt(sqrt(4 * M_PI) * real(Theta2lm(0, 0)));

    const T h = real(hlm(0, 0)) / sqrt(4 * M_PI);

    const T cx = res.cx;
    const T cy = res.cy;
    const T cz = res.cz;
    const T R = sqrt(res.area / (4 * M_PI));

    CCTK_VINFO("iter=%d h=%f c=[%f,%f,%f] R=%f", iter, double(h), double(cx),
               double(cy), double(cz), double(R));
    CCTK_VINFO("  h_0m=%f", double(real(hlm(0, 0))));
    CCTK_VINFO("  h_1m=%f (%f,%f)", double(real(hlm(1, 0))),
               double(real(hlm(1, 1))), double(imag(hlm(1, 1))));
    CCTK_VINFO("  h_2m=%f (%f,%f) (%f,%f)", double(real(hlm(2, 0))),
               double(real(hlm(2, 1))), double(imag(hlm(2, 1))),
               double(real(hlm(2, 2))), double(imag(hlm(2, 2))));
    CCTK_VINFO("  |Θ|∞=%g |Θ|2=%g |Δh|∞=%g |Δh|2=%g", double(Theta_maxabs),
               double(Theta_norm2), double(dh_maxabs), double(dh_norm2));

    const T eps = pow(numeric_limits<T>::epsilon(), T(3) / 4);
    if (iter >= maxiters || dh_maxabs <= eps)
      return res;
    hlm_ptr = make_unique<alm_t<std::complex<T> > >(move(hlm_new));
  }
}

#endif

////////////////////////////////////////////////////////////////////////////////

extern "C" void AHFinder_find(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_find;
  DECLARE_CCTK_PARAMETERS;

  const geom_t geom(npoints);
  if (0)
    std::cout << "geom:" << geom << "\n";

  const auto hlm = scalar_from_const(geom, CCTK_COMPLEX(r0), CCTK_COMPLEX(r1z));
  if (0)
    std::cout << "hlm:" << hlm << "\n";

  const auto hij = evaluate(hlm);
  if (0)
    std::cout << "hij:" << hij << "\n";
  const auto coords = coords_from_shape(hij);
  CCTK_VWARN(CCTK_WARN_ALERT, "Using Brill-Lindquist metric");
  const auto metric = brill_lindquist_metric(cctkGH, coords);
  // const auto metric = interpolate_metric(cctkGH, coords);
  const auto Thetalm = expansion(metric, hlm);
  const auto Thetaij = evaluate(Thetalm);

  const auto Theta_maxabs = maxabs(Thetaij());
  using std::sqrt;
  const auto Theta_avg = real(Thetalm()(0, 0)) / sqrt(4 * M_PI);
  CCTK_VINFO("Θ=%.17g   |Θ|=%.17g", Theta_avg, Theta_maxabs);

#if 0
  const auto res = solve(cctkGH, hlm);
#endif
}

} // namespace AHFinder
