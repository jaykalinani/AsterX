#ifndef DISCRETIZATION_HXX
#define DISCRETIZATION_HXX

#include <cctk.h>

#include <ssht/ssht.h>

#include <array>
#include <cassert>
#include <complex>
#include <vector>

namespace AHFinder {
using namespace std;

struct geom_t {
  const int nmodes;

  // Grid
  const int ntheta;
  const int nphi;
  const int npoints;

  // Coefficients
  const int lmax;
  const int ncoeffs;

  geom_t() = delete;
  geom_t(const int nmodes)
      : nmodes(nmodes), ntheta(nmodes), nphi(2 * nmodes - 1),
        npoints(ntheta * nphi), lmax(nmodes - 1), ncoeffs(nmodes * nmodes) {}

  double coord_theta(const int i, const int j) const {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double theta = ssht_sampling_mw_t2theta(i, nmodes);
    assert(theta >= 0 && theta <= M_PI);
    return theta;
  }

  double coord_phi(const int i, const int j) const {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double phi = ssht_sampling_mw_p2phi(j, nmodes);
    assert(phi >= 0 && phi < 2 * M_PI);
    return phi;
  }

  double coord_dtheta(const int i, const int j) const {
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const double theta_m =
        i == 0 ? 0.0 : (coord_theta(i - 1, j) + coord_theta(i, j)) / 2;
    const double theta_p =
        i == ntheta - 1 ? M_PI
                        : (coord_theta(i, j) + coord_theta(i + 1, j)) / 2;
    return theta_p - theta_m;
  }

  double coord_dphi(const int i, const int j) const { return 2 * M_PI / nphi; }

  int gind(const int i, const int j) const {
    // 0 <= i < ntheta
    // 0 <= j < nphi
    assert(i >= 0 && i < ntheta);
    assert(j >= 0 && j < nphi);
    const int ind = j + nphi * i;
    assert(ind >= 0 && ind < npoints);
    return ind;
  }

  int cind(const int l, const int m) const {
    // 0 <= l <= lmax
    // -l <= m <= l
    assert(l >= 0 && l <= lmax);
    assert(m >= -l && m <= l);
    int ind;
    ssht_sampling_elm2ind(&ind, l, m);
    assert(ind >= 0 && ind < ncoeffs);
    return ind;
  }
};

template <typename T> struct alm_t {
  const geom_t &geom;

  vector<complex<T> > alm;

  alm_t() = delete;
  alm_t(const geom_t &geom) : geom(geom), alm(geom.ncoeffs) {}

  const complex<T> *data() const { return alm.data(); }
  complex<T> *data() { return alm.data(); }
  const complex<T> &operator()(const int l, const int m) const {
    return alm.at(geom.cind(l, m));
  }
  complex<T> &operator()(const int l, const int m) {
    return alm.at(geom.cind(l, m));
  }
};

template <typename T> struct aij_t {
  const geom_t &geom;

  vector<T> aij;

  aij_t() = delete;
  aij_t(const geom_t &geom) : geom(geom), aij(geom.npoints, NAN) {}

  const T *data() const { return aij.data(); }
  T *data() { return aij.data(); }
  const T &operator()(const int i, const int j) const {
    return aij.at(geom.gind(i, j));
  }
  T &operator()(const int i, const int j) { return aij.at(geom.gind(i, j)); }
};

template <typename T>
alm_t<T> coefficients_from_const(const geom_t &geom, const T r0,
                                 const T r1z = 0) {
  alm_t<T> alm(geom);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      alm(l, m) = l == 0 && m == 0 ? sqrt(T(4 * M_PI)) * r0 : 0;
  alm(1, 0) = 2 * r1z;
  return alm;
}

template <typename T> alm_t<T> expand(const aij_t<T> &aij) {
  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]
  alm_t<T> alm(aij.geom);
  ssht_core_mw_forward_sov_conv_sym_real(alm.data(), aij.data(),
                                         alm.geom.nmodes, method, verbosity);
  return alm;
}

template <typename T>
alm_t<T> expand(const aij_t<complex<T> > &aij, const int spin) {
  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]
  alm_t<T> alm(aij.geom);
  ssht_core_mw_forward_sov_conv_sym(alm.data(), aij.data(), alm.geom.nmodes,
                                    spin, method, verbosity);
  return alm;
}

template <typename T> aij_t<T> evaluate(const alm_t<T> &alm) {
  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]
  aij_t<T> aij(alm.geom);
  ssht_core_mw_inverse_sov_sym_real(aij.data(), alm.data(), aij.geom.nmodes,
                                    method, verbosity);
  return aij;
}

template <typename T>
aij_t<complex<T> > evaluate(const alm_t<T> &alm, const int spin) {
  const ssht_dl_method_t method = SSHT_DL_RISBO;
  const int verbosity = 0; // [0..5]
  const geom_t &geom = alm.geom;
  aij_t<complex<T> > aij(geom);
  ssht_core_mw_inverse_sov_sym(aij.data(), alm.data(), aij.geom.nmodes, spin,
                               method, verbosity);
  return aij;
}

template <typename T> alm_t<T> filter(const alm_t<T> &alm, const int lmax) {
  const geom_t &geom = alm.geom;
  alm_t<T> blm(geom);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = l <= lmax ? alm(l, m) : 0;
  return blm;
}

// \dh
template <typename T> alm_t<T> deriv(const alm_t<T> &alm, const int s = 0) {
  const geom_t &geom = alm.geom;
  alm_t<T> blm(geom);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = +sqrt(T((l - s) * (l + s + 1))) * alm(l, m);
  return blm;
}

// \bar\dh
template <typename T> alm_t<T> deriv_bar(const alm_t<T> &alm, const int s = 0) {
  const geom_t &geom = alm.geom;
  alm_t<T> blm(geom);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = -sqrt(T((l + s) * (l - s + 1))) * alm(l, m);
  return blm;
}

template <typename T> alm_t<T> grad(const alm_t<T> &alm) {
  const geom_t &geom = alm.geom;
  alm_t<T> blm(geom);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = -sqrt(T(l * (l + 1))) * alm(l, m);
  return blm;
}

template <typename T> alm_t<T> div(const alm_t<T> &alm) {
  const geom_t &geom = alm.geom;
  alm_t<T> blm(geom);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = sqrt(T(l * (l + 1))) * alm(l, m);
  return blm;
}

// laplace = div grad
template <typename T> alm_t<T> laplace(const alm_t<T> &alm) {
  const geom_t &geom = alm.geom;
  alm_t<T> blm(geom);
  for (int l = 0; l <= geom.lmax; ++l)
    for (int m = -l; m <= l; ++m)
      blm(l, m) = -T(l * (l + 1)) * alm(l, m);
  return blm;
}

template <typename T> alm_t<T> expand_grad(const aij_t<complex<T> > &aij) {
  const geom_t &geom = aij.geom;
  aij_t<complex<T> > aij_p(geom);
  // aij_t<complex<T> > aij_m(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
    for (int j = 0; j < geom.nphi; ++j) {
      aij_p(i, j) = aij(i, j);
      // aij_m(i, j) = conj(aij(i, j));
    }
  }
  const alm_t<T> alm_p = expand(aij_p, +1);
  // const alm_t<T> alm_m = expand(aij_m, -1);
  alm_t<T> alm(geom);
  for (int l = 0; l <= geom.lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      // alm(l, m) = 1 / T(2) * (alm_p(l, m) - alm_m(l, m));
      // const complex<T> blm = -1i / T(2) * (alm_p(l, m) + alm_m(l, m));
      // assert(abs(blm) <= 1.0e-12);
      alm(l, m) = alm_p(l, m);
    }
  }
  return alm;
}

template <typename T> aij_t<complex<T> > evaluate_grad(const alm_t<T> &alm) {
  const geom_t &geom = alm.geom;
  alm_t<T> alm_p(geom);
  // alm_t<T> alm_m(geom);
  for (int l = 0; l <= geom.lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      alm_p(l, m) = alm(l, m);
      // alm_m(l, m) = -alm(l, m);
    }
  }
  const aij_t<complex<T> > aij_p = evaluate(alm_p, +1);
  // const aij_t<complex<T> > aij_m = evaluate(alm_m, -1);
  aij_t<complex<T> > aij(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
    for (int j = 0; j < geom.nphi; ++j) {
      // aij(i, j) = aij_p(i, j) + conj(aij_m(i, j));
      aij(i, j) = aij_p(i, j);
    }
  }
  return aij;
}

} // namespace AHFinder

#endif // #ifndef DISCRETIZATION_HXX
