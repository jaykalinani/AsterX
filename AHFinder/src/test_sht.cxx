#include "sYlm.hxx"

#include <ssht.h>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

namespace AHFinder {
using namespace std;
using namespace std::complex_literals;

extern "C" void AHFinder_test_sht(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_test_sht;

  // Choose resolution

  const int nmodes = 5;

  // Access grid and spectral information

  const int ntheta = nmodes;
  const int nphi = 2 * nmodes - 1;
  const int npoints = ntheta * nphi;
  const auto coord_theta{[&](const int i, const int j) -> CCTK_REAL {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    return ssht_sampling_mw_t2theta(i, nmodes);
  }};
  const auto coord_phi{[&](const int i, const int j) -> CCTK_REAL {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    return ssht_sampling_mw_p2phi(j, nmodes);
  }};
  const auto gind{[&](const int i, const int j) -> int {
    // 0 <= i < ntheta
    // 0 <= j < nphi
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const int ind = j + nphi * i;
    assert(ind >= 0 && ind < npoints);
    return ind;
  }};

  const int lmax = nmodes - 1;
  const int ncoeffs = nmodes * nmodes;
  const auto cind{[&](const int l, const int m) -> int {
    // 0 <= l <= lmax
    // -l <= m <= l
    assert(l >= 0 && l <= lmax);
    assert(m >= -l && m <= l);
    int ind;
    ssht_sampling_elm2ind(&ind, l, m);
    assert(ind >= 0 && ind < ncoeffs);
    return ind;
  }};

  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]

  // Test complex spherical harmonic functions

  for (int spin = sYlm_smin; spin <= sYlm_smax; ++spin) {
    for (int ll = abs(spin); ll <= sYlm_lmax; ++ll) {
      for (int mm = -ll; mm <= ll; ++mm) {

        vector<complex<double> > alm(ncoeffs, NAN);
        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            const complex<double> a = l == ll && m == mm ? 1 : 0;
            alm.at(cind(l, m)) = a;
          }
        }

        vector<complex<double> > aij(npoints, NAN);
        ssht_core_mw_inverse_sov_sym(aij.data(), alm.data(), nmodes, spin,
                                     method, verbosity);

        for (int i = 0; i < ntheta; ++i) {
          for (int j = 0; j < nphi; ++j) {
            const double theta = coord_theta(i, j);
            const double phi = coord_phi(i, j);
            const complex<double> s = sYlm(spin, ll, mm, theta, phi);
            const complex<double> a = aij.at(gind(i, j));
            if (!(fabs(a - s) <= 1.0e-12))
              CCTK_VINFO("spin=%d,ll=%d,mm=%d i=%d,j=%d,theta=%f,phi=%f "
                         "s=(%.17g,%.17g) a=(%.17g,%.17g)",
                         spin, ll, mm, i, j, theta, phi, real(s), imag(s),
                         real(a), imag(a));
            assert(fabs(a - s) <= 1.0e-12);
          }
        }

        ssht_core_mw_forward_sov_conv_sym(alm.data(), aij.data(), nmodes, spin,
                                          method, verbosity);

        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            const complex<double> a = l == ll && m == mm ? 1 : 0;
            assert(fabs(alm.at(cind(l, m)) - a) <= 1.0e-12);
          }
        }
      }
    }
  }

  // Test real spherical harmonic functions

  {
    const int spin = 0;
    for (int ll = abs(spin); ll <= sYlm_lmax; ++ll) {
      for (int mm = -ll; mm <= ll; ++mm) {

        vector<complex<double> > alm(ncoeffs, complex<double>{NAN, NAN});
        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            // f*_lm = (-1)^m f_l,-m
            const complex<double> c = mm < 0 ? 1i : 1;
            complex<double> a{0};
            if (l == ll && m == mm)
              a += c;
            if (l == ll && -m == mm)
              a += double(bitsign(mm)) * conj(c);
            alm.at(cind(l, m)) = a;
          }
        }

        vector<double> aij(npoints, NAN);
        ssht_core_mw_inverse_sov_sym_real(aij.data(), alm.data(), nmodes,
                                          method, verbosity);

        for (int i = 0; i < ntheta; ++i) {
          for (int j = 0; j < nphi; ++j) {
            const double theta = coord_theta(i, j);
            const double phi = coord_phi(i, j);
            const complex<double> c = mm < 0 ? 1i : 1;
            const complex<double> scp = sYlm(spin, ll, mm, theta, phi);
            const complex<double> scm = sYlm(spin, ll, -mm, theta, phi);
            const complex<double> sc =
                c * scp + double(bitsign(mm)) * conj(c) * scm;
            assert(fabs(imag(sc)) <= 1.0e-12);
            const double s = real(sc);
            const double a = aij.at(gind(i, j));
            if (!(fabs(a - s) <= 1.0e-12))
              CCTK_VINFO("spin=%d,ll=%d,mm=%d i=%d,j=%d,theta=%f,phi=%f "
                         "s=%.17g a=%.17g",
                         spin, ll, mm, i, j, theta, phi, s, a);
            assert(fabs(a - s) <= 1.0e-12);
          }
        }

        fill(alm.begin(), alm.end(), complex<double>{NAN, NAN});
        ssht_core_mw_forward_sov_conv_sym_real(alm.data(), aij.data(), nmodes,
                                               method, verbosity);

        for (int l = abs(spin); l <= lmax; ++l) {
          for (int m = -l; m <= l; ++m) {
            const complex<double> c = mm < 0 ? 1i : 1;
            complex<double> a{0};
            if (l == ll && m == mm)
              a += c;
            if (l == ll && -m == mm)
              a += double(bitsign(mm)) * conj(c);
            assert(fabs(alm.at(cind(l, m)) - a) <= 1.0e-12);
          }
        }
      }
    }
  }

  // A note on derivatives: (see
  // <https://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics>:
  //
  // \Delta = \dh \bar\dh = \bar\dh \dh
  //
  //     \dh sYlm = + \sqrt{ (l-s) (l+s+1) (s+1)Ylm
  // \bar\dh sYlm = - \sqrt{ (l+s) (l-s+1) (s-1)Ylm
}

} // namespace AHFinder
