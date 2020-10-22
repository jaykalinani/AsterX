#include "discretization.hxx"
#include "sYlm.hxx"

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <limits>

namespace AHFinder {

template <typename T> constexpr bool is_approx(const T x, const T y) {
  typedef decltype(abs(declval<T>())) R;
  R eps = pow(numeric_limits<R>::epsilon(), R(3) / R(4));
  return abs(x - y) <= eps;
}

extern "C" void AHFinder_test_discretization(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_test_discretization;

  const int nmodes = 5;
  const geom_t geom(nmodes);

  typedef CCTK_REAL T;

  // Loop over all modes
  int count = 0;
  for (int ll = 0; ll < nmodes; ++ll) {
    for (int mm = 0; mm <= ll; ++mm) {
      for (int cc = 0; cc <= min(mm, 1); ++cc) { // re / im
        ++count;

        // Define coefficients (for a real function)
        alm_t<T> alm(geom);
        for (int l = 0; l <= geom.lmax; ++l)
          for (int m = 0; m <= l; ++m)
            alm(l, m) = 0;
        if (mm == 0) {
          alm(ll, mm) = 1;
        } else {
          alm(ll, mm) = cc == 0 ? 1 : 1i;
          alm(ll, -mm) = T(bitsign(mm)) * conj(alm(ll, mm));
        }

        // Evaluate
        const aij_t<T> aij = evaluate(alm);
        // const aij_t<complex<T> > aij = evaluate(alm, 0);

        // Check
        if (ll <= sYlm_lmax) {
          for (int i = 0; i < geom.ntheta; ++i) {
            for (int j = 0; j < geom.nphi; ++j) {
              const T theta = geom.coord_theta(i, j);
              const T phi = geom.coord_phi(i, j);
              const T a = aij(i, j);
              complex<T> sc;
              if (mm == 0) {
                sc = sYlm(0, ll, mm, theta, phi);
              } else {
                const complex<T> c = cc == 0 ? 1 : 1i;
                sc = c * sYlm(0, ll, mm, theta, phi) +
                     T(bitsign(mm)) * conj(c) * sYlm(0, ll, -mm, theta, phi);
              }
              assert(is_approx(imag(sc), T(0)));
              const T s = real(sc);
              assert(is_approx(a, s));
            }
          }
        }

        // Expand
        const alm_t<T> alm1 = expand(aij);
        // const alm_t<T> alm1 = expand(aij, 0);

        // Check
        for (int l = 0; l <= geom.lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            if (!(is_approx(alm1(l, m), alm(l, m))))
              CCTK_VINFO("ll=%d,mm=%d,cc=%d l=%d,m=%d "
                         "alm=(%.17g,%.17g) alm1=(%.17g,%.17g)",
                         ll, mm, cc, l, m, real(alm(l, m)), imag(alm(l, m)),
                         real(alm1(l, m)), imag(alm1(l, m)));
            assert(is_approx(alm1(l, m), alm(l, m)));
          }
        }

        // Gradient
        const alm_t<T> dalm = grad(alm);

        // Evaluate
        const aij_t<complex<T> > daij = evaluate_grad(dalm);

        // Check
        if (ll <= sYlm_lmax) {
          for (int i = 0; i < geom.ntheta; ++i) {
            for (int j = 0; j < geom.nphi; ++j) {
              const T theta = geom.coord_theta(i, j);
              const T phi = geom.coord_phi(i, j);
              const complex<T> da = daij(i, j);
              array<complex<T>, 2> dsc;
              if (mm == 0) {
                dsc = dsYlm(0, ll, mm, theta, phi);
              } else {
                const complex<T> c = cc == 0 ? 1 : 1i;
                const array<complex<T>, 2> dsc_p = dsYlm(0, ll, mm, theta, phi);
                const array<complex<T>, 2> dsc_m =
                    dsYlm(0, ll, -mm, theta, phi);
                for (int n = 0; n < 2; ++n) {
                  dsc[n] = c * dsc_p[n] + T(bitsign(mm)) * conj(c) * dsc_m[n];
                  assert(is_approx(imag(dsc[n]), T(0)));
                }
              }
              complex<T> ds{real(dsc[0]), real(dsc[1])};
              if (!(is_approx(da, ds)))
                CCTK_VINFO("ll=%d,mm=%d,cc=%d i=%d,j=%d,theta=%f,phi=%f "
                           "da=(%f,%f) ds=(%f,%f)",
                           ll, mm, cc, i, j, theta, phi, real(da), imag(da),
                           real(ds), imag(ds));
              assert(is_approx(da, ds));
            }
          }
        }

        // Expand
        const alm_t<T> dalm1 = expand_grad(daij);

        // Check
        for (int l = 0; l <= geom.lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            if (!(is_approx(dalm1(l, m), dalm(l, m))))
              CCTK_VINFO("ll=%d,mm=%d,cc=%d l=%d,m=%d "
                         "dalm=(%.17g,%.17g) dalm1=(%.17g,%.17g)",
                         ll, mm, cc, l, m, real(dalm(l, m)), imag(dalm(l, m)),
                         real(dalm1(l, m)), imag(dalm1(l, m)));
            assert(is_approx(dalm1(l, m), dalm(l, m)));
          }
        }

        // Divergence
        const alm_t<T> lalm = div(dalm);

        // Check
        for (int l = 0; l <= geom.lmax; ++l)
          for (int m = -l; m <= l; ++m)
            assert(is_approx(lalm(l, m), T(-l * (l + 1)) * alm(l, m)));
      }
    }
  }
  assert(count == nmodes * nmodes);
}

} // namespace AHFinder
