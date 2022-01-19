#include "discretization.hxx"

#include <dual.hxx>
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
using namespace std;
using namespace Arith;

template <typename T> struct coords_t {
  geom_t geom;
  aij_t<T> x, y, z;

  coords_t() = delete;
  coords_t(const geom_t &geom) : geom(geom), x(geom), y(geom), z(geom) {}
};

template <typename T> coords_t<T> coords_from_shape(const aij_t<T> &h) {
  DECLARE_CCTK_PARAMETERS;
  const geom_t &geom = h.geom;
  coords_t<T> coords(geom);
  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {
      const T r = h(i, j);
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);
      coords.x(i, j) = x0 + r * sin(theta) * cos(phi);
      coords.y(i, j) = y0 + r * sin(theta) * sin(phi);
      coords.z(i, j) = z0 + r * cos(theta);
    }
  }
  return coords;
}

template <typename T> struct metric_t {
  geom_t geom;
  aij_t<T> gxx, gxy, gxz, gyy, gyz, gzz;
  aij_t<T> gxx_x, gxy_x, gxz_x, gyy_x, gyz_x, gzz_x;
  aij_t<T> gxx_y, gxy_y, gxz_y, gyy_y, gyz_y, gzz_y;
  aij_t<T> gxx_z, gxy_z, gxz_z, gyy_z, gyz_z, gzz_z;
  aij_t<T> kxx, kxy, kxz, kyy, kyz, kzz;

  metric_t() = delete;
  metric_t(const geom_t &geom)
      : geom(geom), gxx(geom), gxy(geom), gxz(geom), gyy(geom), gyz(geom),
        gzz(geom), gxx_x(geom), gxy_x(geom), gxz_x(geom), gyy_x(geom),
        gyz_x(geom), gzz_x(geom), gxx_y(geom), gxy_y(geom), gxz_y(geom),
        gyy_y(geom), gyz_y(geom), gzz_y(geom), gxx_z(geom), gxy_z(geom),
        gxz_z(geom), gyy_z(geom), gyz_z(geom), gzz_z(geom), kxx(geom),
        kxy(geom), kxz(geom), kyy(geom), kyz(geom), kzz(geom) {}
};

template <typename T>
metric_t<T> brill_lindquist_metric(const cGH *const cctkGH,
                                   const coords_t<T> &coords) {
  const geom_t &geom = coords.geom;
  metric_t<T> metric(geom);

  const T x0 = 0.0;
  const T y0 = 0.0;
  const T z0 = 0.0;
  const T fx = 0.8;
  const T fy = 1.0;
  const T fz = 1.0;
  const T M = 1.0;

#pragma omp simd
  for (int n = 0; n < geom.npoints; ++n) {
    using U = vect<T, 3>;
    using DT = dual<T, U>; // spatial derivatives (x, y, z)
    const DT x{coords.x.data()[n], {1, 0, 0}};
    const DT y{coords.y.data()[n], {0, 1, 0}};
    const DT z{coords.z.data()[n], {0, 0, 1}};
    const DT X = (x - x0) / fx;
    const DT Y = (y - y0) / fy;
    const DT Z = (z - z0) / fz;
    const DT R = sqrt(pow(X, 2) + pow(Y, 2) + pow(Z, 2));
    const DT Psi = pow(1 + M / (2 * R), 4);
    const vect<T, 3> eta{pow(X.eps[0], 2), pow(Y.eps[1], 2), pow(Z.eps[2], 2)};
    metric.gxx.data()[n] = Psi.val * eta[0];
    metric.gxy.data()[n] = 0;
    metric.gxz.data()[n] = 0;
    metric.gyy.data()[n] = Psi.val * eta[1];
    metric.gyz.data()[n] = 0;
    metric.gzz.data()[n] = Psi.val * eta[2];
    metric.gxx_x.data()[n] = Psi.eps[0] * eta[0];
    metric.gxy_x.data()[n] = 0;
    metric.gxz_x.data()[n] = 0;
    metric.gyy_x.data()[n] = Psi.eps[0] * eta[1];
    metric.gyz_x.data()[n] = 0;
    metric.gzz_x.data()[n] = Psi.eps[0] * eta[2];
    metric.gxx_y.data()[n] = Psi.eps[1] * eta[0];
    metric.gxy_y.data()[n] = 0;
    metric.gxz_y.data()[n] = 0;
    metric.gyy_y.data()[n] = Psi.eps[1] * eta[1];
    metric.gyz_y.data()[n] = 0;
    metric.gzz_y.data()[n] = Psi.eps[1] * eta[2];
    metric.gxx_z.data()[n] = Psi.eps[2] * eta[0];
    metric.gxy_z.data()[n] = 0;
    metric.gxz_z.data()[n] = 0;
    metric.gyy_z.data()[n] = Psi.eps[2] * eta[1];
    metric.gyz_z.data()[n] = 0;
    metric.gzz_z.data()[n] = Psi.eps[2] * eta[2];
    metric.kxx.data()[n] = 0;
    metric.kxy.data()[n] = 0;
    metric.kxz.data()[n] = 0;
    metric.kyy.data()[n] = 0;
    metric.kyz.data()[n] = 0;
    metric.kzz.data()[n] = 0;
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
      metric.gxx.data(),   metric.gxy.data(),   metric.gxz.data(),
      metric.gyy.data(),   metric.gyz.data(),   metric.gzz.data(),
      metric.gxx_x.data(), metric.gxy_x.data(), metric.gxz_x.data(),
      metric.gyy_x.data(), metric.gyz_x.data(), metric.gzz_x.data(),
      metric.gxx_y.data(), metric.gxy_y.data(), metric.gxz_y.data(),
      metric.gyy_y.data(), metric.gyz_y.data(), metric.gzz_y.data(),
      metric.gxx_z.data(), metric.gxy_z.data(), metric.gxz_z.data(),
      metric.gyy_z.data(), metric.gyz_z.data(), metric.gzz_z.data(),
      metric.kxx.data(),   metric.kxy.data(),   metric.kxz.data(),
      metric.kyy.data(),   metric.kyz.data(),   metric.kzz.data()};

  Interpolate(cctkGH, geom.npoints, coords.x.data(), coords.y.data(),
              coords.z.data(), nvars, varinds.data(), operations.data(),
              ptrs.data());

  return metric;
}

////////////////////////////////////////////////////////////////////////////////

// Expansion and updated horizon shape; see [arXiv:gr-qc/0702038], (29)
template <typename T> struct expansion_t {
  alm_t<T> hlm;
  T area;
  T cx, cy, cz;
  alm_t<T> Thetalm, hlm_new;
};

template <typename T>
expansion_t<T> expansion(const metric_t<T> &metric, const alm_t<T> &hlm) {
  DECLARE_CCTK_PARAMETERS;

  const geom_t &geom = hlm.geom;

  // Evaluate h^ij and its derivatives
  const aij_t<T> hij = evaluate(hlm);
  const alm_t<T> dhlm = grad(hlm);
  const aij_t<complex<T> > dhij = evaluate_grad(dhlm);

  aij_t<T> surij(geom), sutij(geom), z_supij(geom);
  aij_t<T> s_dsur_rij(geom);
  aij_t<complex<T> > s_dsuij(geom);

  aij_t<T> lambdaij(geom);

  aij_t<T> cxij(geom), cyij(geom), czij(geom);
  aij_t<T> Aij(geom);

  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {

      // Coordinates

      const T h = hij(i, j);
      // dX = d/d\theta X + i/\sin\theta d/d\phi X
      const complex<T> dh = dhij(i, j);
      const T h_t = real(dh);
      const T s_h_p = imag(dh);

      // Dual quantities for radial derivatives
      const dual<T> r{h, 1};
      const dual<T> theta{geom.coord_theta(i, j), 0};
      const dual<T> phi{geom.coord_phi(i, j), 0};

      const dual<T> x = x0 + r * sin(theta) * cos(phi);
      const dual<T> y = y0 + r * sin(theta) * sin(phi);
      const dual<T> z = z0 + r * cos(theta);

      const dual<T> x_r = sin(theta) * cos(phi);
      const dual<T> y_r = sin(theta) * sin(phi);
      const dual<T> z_r = cos(theta);
      const dual<T> x_t = r * cos(theta) * cos(phi);
      const dual<T> y_t = r * cos(theta) * sin(phi);
      const dual<T> z_t = -r * sin(theta);
      const dual<T> s_x_p = -r * sin(phi);
      const dual<T> s_y_p = r * cos(phi);
      const dual<T> s_z_p = 0;

      // Metric

      const T gxx0 = metric.gxx(i, j);
      const T gxy0 = metric.gxy(i, j);
      const T gxz0 = metric.gxz(i, j);
      const T gyy0 = metric.gyy(i, j);
      const T gyz0 = metric.gyz(i, j);
      const T gzz0 = metric.gzz(i, j);

      const T gxx0_x = metric.gxx_x(i, j);
      const T gxy0_x = metric.gxy_x(i, j);
      const T gxz0_x = metric.gxz_x(i, j);
      const T gyy0_x = metric.gyy_x(i, j);
      const T gyz0_x = metric.gyz_x(i, j);
      const T gzz0_x = metric.gzz_x(i, j);
      const T gxx0_y = metric.gxx_y(i, j);
      const T gxy0_y = metric.gxy_y(i, j);
      const T gxz0_y = metric.gxz_y(i, j);
      const T gyy0_y = metric.gyy_y(i, j);
      const T gyz0_y = metric.gyz_y(i, j);
      const T gzz0_y = metric.gzz_y(i, j);
      const T gxx0_z = metric.gxx_z(i, j);
      const T gxy0_z = metric.gxy_z(i, j);
      const T gxz0_z = metric.gxz_z(i, j);
      const T gyy0_z = metric.gyy_z(i, j);
      const T gyz0_z = metric.gyz_z(i, j);
      const T gzz0_z = metric.gzz_z(i, j);

      // Radial derivative of metric
      const T gxx0_r = gxx0_x * x_r.val + gxx0_y * y_r.val + gxx0_z * z_r.val;
      const T gxy0_r = gxy0_x * x_r.val + gxy0_y * y_r.val + gxy0_z * z_r.val;
      const T gxz0_r = gxz0_x * x_r.val + gxz0_y * y_r.val + gxz0_z * z_r.val;
      const T gyy0_r = gyy0_x * x_r.val + gyy0_y * y_r.val + gyy0_z * z_r.val;
      const T gyz0_r = gyz0_x * x_r.val + gyz0_y * y_r.val + gyz0_z * z_r.val;
      const T gzz0_r = gzz0_x * x_r.val + gzz0_y * y_r.val + gzz0_z * z_r.val;

      const dual<T> gxx{gxx0, gxx0_r};
      const dual<T> gxy{gxy0, gxy0_r};
      const dual<T> gxz{gxz0, gxz0_r};
      const dual<T> gyy{gyy0, gyy0_r};
      const dual<T> gyz{gyz0, gyz0_r};
      const dual<T> gzz{gzz0, gzz0_r};

      // Metric in spherical coordinates
      const dual<T> qrr = 0                 //
                          + gxx * x_r * x_r //
                          + gxy * x_r * y_r //
                          + gxz * x_r * z_r //
                          + gxy * y_r * x_r //
                          + gyy * y_r * y_r //
                          + gyz * y_r * z_r //
                          + gxz * z_r * x_r //
                          + gyz * z_r * y_r //
                          + gzz * z_r * z_r;
      const dual<T> qrt = 0                 //
                          + gxx * x_r * x_t //
                          + gxy * x_r * y_t //
                          + gxz * x_r * z_t //
                          + gxy * y_r * x_t //
                          + gyy * y_r * y_t //
                          + gyz * y_r * z_t //
                          + gxz * z_r * x_t //
                          + gyz * z_r * y_t //
                          + gzz * z_r * z_t;
      const dual<T> s_qrp = 0                   //
                            + gxx * x_r * s_x_p //
                            + gxy * x_r * s_y_p //
                            + gxz * x_r * s_z_p //
                            + gxy * y_r * s_x_p //
                            + gyy * y_r * s_y_p //
                            + gyz * y_r * s_z_p //
                            + gxz * z_r * s_x_p //
                            + gyz * z_r * s_y_p //
                            + gzz * z_r * s_z_p;
      const dual<T> qtt = 0                 //
                          + gxx * x_t * x_t //
                          + gxy * x_t * y_t //
                          + gxz * x_t * z_t //
                          + gxy * y_t * x_t //
                          + gyy * y_t * y_t //
                          + gyz * y_t * z_t //
                          + gxz * z_t * x_t //
                          + gyz * z_t * y_t //
                          + gzz * z_t * z_t;
      const dual<T> s_qtp = 0                   //
                            + gxx * x_t * s_x_p //
                            + gxy * x_t * s_y_p //
                            + gxz * x_t * s_z_p //
                            + gxy * y_t * s_x_p //
                            + gyy * y_t * s_y_p //
                            + gyz * y_t * s_z_p //
                            + gxz * z_t * s_x_p //
                            + gyz * z_t * s_y_p //
                            + gzz * z_t * s_z_p;
      const dual<T> ss_qpp = 0                     //
                             + gxx * s_x_p * s_x_p //
                             + gxy * s_x_p * s_y_p //
                             + gxz * s_x_p * s_z_p //
                             + gxy * s_y_p * s_x_p //
                             + gyy * s_y_p * s_y_p //
                             + gyz * s_y_p * s_z_p //
                             + gxz * s_z_p * s_x_p //
                             + gyz * s_z_p * s_y_p //
                             + gzz * s_z_p * s_z_p;

      const dual<T> ss_detq = 0                                      //
                              + qrr * (qtt * ss_qpp - pow(s_qtp, 2)) //
                              + qrt * (s_qtp * s_qrp - qrt * ss_qpp) //
                              + s_qrp * (qrt * s_qtp - qtt * s_qrp);
      const dual<T> s_sqrt_detq = sqrt(ss_detq);

      const dual<T> qurr = (qtt * ss_qpp - pow(s_qtp, 2)) / ss_detq;
      const dual<T> qurt = (s_qtp * s_qrp - ss_qpp * qrt) / ss_detq;
      const dual<T> z_qurp = (qrt * s_qtp - s_qrp * qtt) / ss_detq;
      const dual<T> qutt = (ss_qpp * qrr - pow(s_qrp, 2)) / ss_detq;
      const dual<T> z_qutp = (s_qrp * qrt - qrr * s_qtp) / ss_detq;
      const dual<T> zz_qupp = (qrr * qtt - pow(qrt, 2)) / ss_detq;

#ifdef CCTK_DEBUG
      const dual<T> q[3][3]{
          {qrr, qrt, s_qrp}, {qrt, qtt, s_qtp}, {s_qrp, s_qtp, ss_qpp}};
      const dual<T> qu[3][3]{{qurr, qurt, z_qurp},
                             {qurt, qutt, z_qutp},
                             {z_qurp, z_qutp, zz_qupp}};
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          assert(q[a][b] == q[b][a]);
          assert(qu[a][b] == qu[b][a]);
          dual<T> s = 0;
          for (int c = 0; c < 3; ++c)
            s += qu[a][c] * q[c][b];
          assert(fabs(s - (a == b)) <= 1.0e-12);
        }
      }
#endif

      // Spacelike normal

      // const T F = r.val - h;

      const dual<T> F_r = {1, 0};
      const dual<T> F_t = {-h_t, 0};
      const dual<T> s_F_p = {-s_h_p, 0};

      const dual<T> Fu_r = qurr * F_r + qurt * F_t + z_qurp * s_F_p;
      const dual<T> Fu_t = qurt * F_r + qutt * F_t + z_qutp * s_F_p;
      const dual<T> z_Fu_p = z_qurp * F_r + z_qutp * F_t + zz_qupp * s_F_p;

      const dual<T> dF2 = F_r * Fu_r + F_t * Fu_t + s_F_p * z_Fu_p;
      const dual<T> dF = sqrt(dF2);

      // spacelike normal s_i
      const dual<T> sr = F_r / dF;
      const dual<T> st = F_t / dF;
      const dual<T> s_sp = s_F_p / dF;

      const dual<T> sur = qurr * sr + qurt * st + z_qurp * s_sp;
      const dual<T> sut = qurt * sr + qutt * st + z_qutp * s_sp;
      const dual<T> z_sup = z_qurp * sr + z_qutp * st + zz_qupp * s_sp;

#ifdef CCTK_DEBUG
      const dual<T> s2 = sr * sur + st * sut + s_sp * z_sup;
      assert(fabs(s2 - 1) <= 1.0e-12);
#endif

#ifdef CCTK_DEBUG
      // m1_a
      // TODO: rename to m1^a
      dual<T> m1ur = h_t;
      dual<T> m1ut = 1;
      dual<T> z_m1up = 0;

      dual<T> m1r = qrr * m1ur + qrt * m1ut + s_qrp * z_m1up;
      dual<T> m1t = qrt * m1ur + qtt * m1ut + s_qtp * z_m1up;
      dual<T> s_m1p = s_qrp * m1ur + s_qtp * m1ut + ss_qpp * z_m1up;

      dual<T> m12 = m1r * m1ur + m1t * m1ut + s_m1p * z_m1up;

      m1ur /= sqrt(m12);
      m1ut /= sqrt(m12);
      z_m1up /= sqrt(m12);

      m1r = qrr * m1ur + qrt * m1ut + s_qrp * z_m1up;
      m1t = qrt * m1ur + qtt * m1ut + s_qtp * z_m1up;
      s_m1p = s_qrp * m1ur + s_qtp * m1ut + ss_qpp * z_m1up;

      m12 = m1r * m1ur + m1t * m1ut + s_m1p * z_m1up;
      assert(fabs(m12 - 1) <= 1.0e-12);

      const dual<T> sm1 = sr * m1ur + st * m1ut + s_sp * z_m1up;
      assert(fabs(sm1) <= 1.0e-12);

      // m2_a
      dual<T> s_m2ur = s_h_p;
      dual<T> s_m2ut = 0;
      dual<T> m2up = 1;

      const dual<T> s_m1m2 = s_m2ur * m1r + s_m2ut * m1t + m2up * s_m1p;

      s_m2ur -= s_m1m2 * m1ur;
      s_m2ut -= s_m1m2 * m1ut;
      m2up -= s_m1m2 * z_m1up;

      const dual<T> s_m2r = qrr * s_m2ur + qrt * s_m2ut + s_qrp * m2up;
      const dual<T> s_m2t = qrt * s_m2ur + qtt * s_m2ut + s_qtp * m2up;
      const dual<T> ss_m2p = s_qrp * s_m2ur + s_qtp * s_m2ut + ss_qpp * m2up;

      const dual<T> ss_m22 = s_m2r * s_m2ur + s_m2t * s_m2ut + ss_m2p * m2up;

      const dual<T> m2ur = s_m2ur / sqrt(ss_m22);
      const dual<T> m2ut = s_m2ut / sqrt(ss_m22);
      const dual<T> z_m2up = m2up / sqrt(ss_m22);

      const dual<T> m2r = qrr * m2ur + qrt * m2ut + s_qrp * z_m2up;
      const dual<T> m2t = qrt * m2ur + qtt * m2ut + s_qtp * z_m2up;
      const dual<T> s_m2p = s_qrp * m2ur + s_qtp * m2ut + ss_qpp * z_m2up;

      const dual<T> m1m2 = m2ur * m1r + m2ut * m1t + z_m2up * s_m1p;
      assert(fabs(m1m2) <= 1.0e-12);

      const dual<T> m22 = m2r * m2ur + m2t * m2ut + s_m2p * z_m2up;
      assert(fabs(m22 - 1) <= 1.0e-12);

      const dual<T> sm2 = sr * m2ur + st * m2ut + s_sp * z_m2up;
      assert(fabs(sm2) <= 1.0e-12);
#endif

#if 0
      {
        // m1_a
        const dual<T> m1ux = x_t + h_t * sin(theta) * cos(phi);
        const dual<T> m1uy = y_t + h_t * sin(theta) * sin(phi);
        const dual<T> m1uz = z_t + h_t * cos(theta);

        const dual<T> m1x = gxx * m1ux + gxy * m1uy + gxz * m1uz;
        const dual<T> m1y = gxy * m1ux + gyy * m1uy + gyz * m1uz;
        const dual<T> m1z = gxz * m1ux + gyz * m1uy + gzz * m1uz;

        const dual<T> m1r = m1x * x_r + m1y * y_r + m1z * z_r;
        const dual<T> m1t = m1x * x_t + m1y * y_t + m1z * z_t;
        const dual<T> s_m1p = m1x * s_x_p + m1y * s_y_p + m1z * s_z_p;

        const dual<T> sm1 = sur * m1r + sut * m1t + z_sup * s_m1p;
        assert(fabs(sm1) <= 1.0e-12);

        // m2_a
        const dual<T> s_m2ux = s_x_p + s_h_p * sin(theta) * cos(phi);
        const dual<T> s_m2uy = s_y_p + s_h_p * sin(theta) * sin(phi);
        const dual<T> s_m2uz = s_z_p + s_h_p * cos(theta);

        const dual<T> s_m2x = gxx * s_m2ux + gxy * s_m2uy + gxz * s_m2uz;
        const dual<T> s_m2y = gxy * s_m2ux + gyy * s_m2uy + gyz * s_m2uz;
        const dual<T> s_m2z = gxz * s_m2ux + gyz * s_m2uy + gzz * s_m2uz;

        const dual<T> s_m2r = s_m2x * x_r + s_m2y * y_r + s_m2z * z_r;
        const dual<T> s_m2t = s_m2x * x_t + s_m2y * y_t + s_m2z * z_t;
        const dual<T> ss_m2p = s_m2x * s_x_p + s_m2y * s_y_p + s_m2z * s_z_p;

        const dual<T> s_sm2 = sur * s_m2r + sut * s_m2t + z_sup * ss_m2p;
        assert(fabs(s_sm2) <= 1.0e-12);
      }
#endif

      // densitized spacelike normal
      const dual<T> s_dsur = s_sqrt_detq * sur;
      const dual<T> s_dsut = s_sqrt_detq * sut;
      const dual<T> dsup = s_sqrt_detq * z_sup;

      const T Psi4 = cbrt(ss_detq.val / pow(r.val, 4));
      // [arXiv:gr-qc/0702038], (26)
      const T lambda = Psi4 * dF.val * pow(r.val, 2);

      // Construct a 2-metric on the surface

      const dual<T> x_T =
          h_t * sin(theta) * cos(phi) + r * cos(theta) * cos(phi);
      const dual<T> y_T =
          h_t * sin(theta) * sin(phi) + r * cos(theta) * sin(phi);
      const dual<T> z_T = h_t * cos(theta) - r * sin(theta);
      const dual<T> s_x_P = s_h_p * sin(theta) * cos(phi) - r * sin(phi);
      const dual<T> s_y_P = s_h_p * sin(theta) * sin(phi) + r * cos(phi);
      const dual<T> s_z_P = s_h_p * cos(theta);

      const dual<T> qTT = 0                 //
                          + gxx * x_T * x_T //
                          + gxy * x_T * y_T //
                          + gxz * x_T * z_T //
                          + gxy * y_T * x_T //
                          + gyy * y_T * y_T //
                          + gyz * y_T * z_T //
                          + gxz * z_T * x_T //
                          + gyz * z_T * y_T //
                          + gzz * z_T * z_T;
      const dual<T> s_qTP = 0                   //
                            + gxx * x_T * s_x_P //
                            + gxy * x_T * s_y_P //
                            + gxz * x_T * s_z_P //
                            + gxy * y_T * s_x_P //
                            + gyy * y_T * s_y_P //
                            + gyz * y_T * s_z_P //
                            + gxz * z_T * s_x_P //
                            + gyz * z_T * s_y_P //
                            + gzz * z_T * s_z_P;
      const dual<T> ss_qPP = 0                     //
                             + gxx * s_x_P * s_x_P //
                             + gxy * s_x_P * s_y_P //
                             + gxz * s_x_P * s_z_P //
                             + gxy * s_y_P * s_x_P //
                             + gyy * s_y_P * s_y_P //
                             + gyz * s_y_P * s_z_P //
                             + gxz * s_z_P * s_x_P //
                             + gyz * s_z_P * s_y_P //
                             + gzz * s_z_P * s_z_P;

      const dual<T> ss_detq2 = qTT * ss_qPP - pow(s_qTP, 2);
      const dual<T> s_sqrt_detq2 = sqrt(ss_detq2);
      const T A = s_sqrt_detq2.val;

      // spacelike normal
      surij(i, j) = sur.val;
      sutij(i, j) = sut.val;
      z_supij(i, j) = z_sup.val;

      // densitized spacelike normal
      // r derivative of r component
      s_dsur_rij(i, j) = s_dsur.eps;
      // theta and phi components (will calculate 2-divergence below)
      s_dsuij(i, j) = complex<T>(s_dsut.val, dsup.val);

      lambdaij(i, j) = lambda;

      Aij(i, j) = A;
      cxij(i, j) = A * x.val;
      cyij(i, j) = A * y.val;
      czij(i, j) = A * z.val;
    }
  }

  const alm_t<T> s_dsulm = expand_grad(s_dsuij);
  const alm_t<T> s_lsulm = div(s_dsulm);
  const aij_t<T> s_lsuij = evaluate(s_lsulm);

  const alm_t<T> Alm = expand(Aij);
  const alm_t<T> cxlm = expand(cxij);
  const alm_t<T> cylm = expand(cyij);
  const alm_t<T> czlm = expand(czij);

  aij_t<T> Thetaij(geom);

  for (int i = 0; i < geom.ntheta; ++i) {
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j) {

      // Coordinates

      const T h = hij(i, j);
      // dX = d/d\theta X + i/\sin\theta d/d\phi X
      // const complex<T> dh = dhij(i, j);
      // const T h_t = real(dh);
      // const T s_h_p = imag(dh);

      const T r = h;
      const T theta = geom.coord_theta(i, j);
      const T phi = geom.coord_phi(i, j);

      // const T x = x0 + r * sin(theta) * cos(phi);
      // const T y = y0 + r * sin(theta) * sin(phi);
      // const T z = z0 + r * cos(theta);

      const T x_r = sin(theta) * cos(phi);
      const T y_r = sin(theta) * sin(phi);
      const T z_r = cos(theta);
      const T x_t = r * cos(theta) * cos(phi);
      const T y_t = r * cos(theta) * sin(phi);
      const T z_t = -r * sin(theta);
      const T s_x_p = -r * sin(phi);
      const T s_y_p = r * cos(phi);
      const T s_z_p = 0;

      // Metric

      const T gxx = metric.gxx(i, j);
      const T gxy = metric.gxy(i, j);
      const T gxz = metric.gxz(i, j);
      const T gyy = metric.gyy(i, j);
      const T gyz = metric.gyz(i, j);
      const T gzz = metric.gzz(i, j);

      // Metric in spherical coordinates
      const T qrr = 0                 //
                    + gxx * x_r * x_r //
                    + gxy * x_r * y_r //
                    + gxz * x_r * z_r //
                    + gxy * y_r * x_r //
                    + gyy * y_r * y_r //
                    + gyz * y_r * z_r //
                    + gxz * z_r * x_r //
                    + gyz * z_r * y_r //
                    + gzz * z_r * z_r;
      const T qrt = 0                 //
                    + gxx * x_r * x_t //
                    + gxy * x_r * y_t //
                    + gxz * x_r * z_t //
                    + gxy * y_r * x_t //
                    + gyy * y_r * y_t //
                    + gyz * y_r * z_t //
                    + gxz * z_r * x_t //
                    + gyz * z_r * y_t //
                    + gzz * z_r * z_t;
      const T s_qrp = 0                   //
                      + gxx * x_r * s_x_p //
                      + gxy * x_r * s_y_p //
                      + gxz * x_r * s_z_p //
                      + gxy * y_r * s_x_p //
                      + gyy * y_r * s_y_p //
                      + gyz * y_r * s_z_p //
                      + gxz * z_r * s_x_p //
                      + gyz * z_r * s_y_p //
                      + gzz * z_r * s_z_p;
      const T qtt = 0                 //
                    + gxx * x_t * x_t //
                    + gxy * x_t * y_t //
                    + gxz * x_t * z_t //
                    + gxy * y_t * x_t //
                    + gyy * y_t * y_t //
                    + gyz * y_t * z_t //
                    + gxz * z_t * x_t //
                    + gyz * z_t * y_t //
                    + gzz * z_t * z_t;
      const T s_qtp = 0                   //
                      + gxx * x_t * s_x_p //
                      + gxy * x_t * s_y_p //
                      + gxz * x_t * s_z_p //
                      + gxy * y_t * s_x_p //
                      + gyy * y_t * s_y_p //
                      + gyz * y_t * s_z_p //
                      + gxz * z_t * s_x_p //
                      + gyz * z_t * s_y_p //
                      + gzz * z_t * s_z_p;
      const T ss_qpp = 0                     //
                       + gxx * s_x_p * s_x_p //
                       + gxy * s_x_p * s_y_p //
                       + gxz * s_x_p * s_z_p //
                       + gxy * s_y_p * s_x_p //
                       + gyy * s_y_p * s_y_p //
                       + gyz * s_y_p * s_z_p //
                       + gxz * s_z_p * s_x_p //
                       + gyz * s_z_p * s_y_p //
                       + gzz * s_z_p * s_z_p;

      const T ss_detq = 0                                      //
                        + qrr * (qtt * ss_qpp - pow(s_qtp, 2)) //
                        + qrt * (s_qtp * s_qrp - qrt * ss_qpp) //
                        + s_qrp * (qrt * s_qtp - qtt * s_qrp);
      const T s_sqrt_detq = sqrt(ss_detq);

      const T qurr = (qtt * ss_qpp - pow(s_qtp, 2)) / ss_detq;
      const T qurt = (s_qtp * s_qrp - ss_qpp * qrt) / ss_detq;
      const T z_qurp = (qrt * s_qtp - s_qrp * qtt) / ss_detq;
      const T qutt = (ss_qpp * qrr - pow(s_qrp, 2)) / ss_detq;
      const T z_qutp = (s_qrp * qrt - qrr * s_qtp) / ss_detq;
      const T zz_qupp = (qrr * qtt - pow(qrt, 2)) / ss_detq;

#ifdef CCTK_DEBUG
      const dual<T> q[3][3]{
          {qrr, qrt, s_qrp}, {qrt, qtt, s_qtp}, {s_qrp, s_qtp, ss_qpp}};
      const dual<T> qu[3][3]{{qurr, qurt, z_qurp},
                             {qurt, qutt, z_qutp},
                             {z_qurp, z_qutp, zz_qupp}};
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) {
          assert(q[a][b] == q[b][a]);
          assert(qu[a][b] == qu[b][a]);
          dual<T> s = 0;
          for (int c = 0; c < 3; ++c)
            s += qu[a][c] * q[c][b];
          assert(fabs(s - (a == b)) <= 1.0e-12);
        }
      }
#endif

      // Spacelike normal

      const T sur = surij(i, j);
      const T sut = sutij(i, j);
      const T z_sup = z_supij(i, j);

      const T s_dsur_r = s_dsur_rij(i, j);
      const T s_lsu = s_lsuij(i, j);

      // Extrinsic curvature

      const T kxx = metric.kxx(i, j);
      const T kxy = metric.kxy(i, j);
      const T kxz = metric.kxz(i, j);
      const T kyy = metric.kyy(i, j);
      const T kyz = metric.kyz(i, j);
      const T kzz = metric.kzz(i, j);

      const T krr = 0                 //
                    + kxx * x_r * x_r //
                    + kxy * x_r * y_r //
                    + kxz * x_r * z_r //
                    + kxy * y_r * x_r //
                    + kyy * y_r * y_r //
                    + kyz * y_r * z_r //
                    + kxz * z_r * x_r //
                    + kyz * z_r * y_r //
                    + kzz * z_r * z_r;
      const T krt = 0                 //
                    + kxx * x_r * x_t //
                    + kxy * x_r * y_t //
                    + kxz * x_r * z_t //
                    + kxy * y_r * x_t //
                    + kyy * y_r * y_t //
                    + kyz * y_r * z_t //
                    + kxz * z_r * x_t //
                    + kyz * z_r * y_t //
                    + kzz * z_r * z_t;
      const T s_krp = 0                   //
                      + kxx * x_r * s_x_p //
                      + kxy * x_r * s_y_p //
                      + kxz * x_r * s_z_p //
                      + kxy * y_r * s_x_p //
                      + kyy * y_r * s_y_p //
                      + kyz * y_r * s_z_p //
                      + kxz * z_r * s_x_p //
                      + kyz * z_r * s_y_p //
                      + kzz * z_r * s_z_p;
      const T ktt = 0                 //
                    + kxx * x_t * x_t //
                    + kxy * x_t * y_t //
                    + kxz * x_t * z_t //
                    + kxy * y_t * x_t //
                    + kyy * y_t * y_t //
                    + kyz * y_t * z_t //
                    + kxz * z_t * x_t //
                    + kyz * z_t * y_t //
                    + kzz * z_t * z_t;
      const T s_ktp = 0                   //
                      + kxx * x_t * s_x_p //
                      + kxy * x_t * s_y_p //
                      + kxz * x_t * s_z_p //
                      + kxy * y_t * s_x_p //
                      + kyy * y_t * s_y_p //
                      + kyz * y_t * s_z_p //
                      + kxz * z_t * s_x_p //
                      + kyz * z_t * s_y_p //
                      + kzz * z_t * s_z_p;
      const T ss_kpp = 0                     //
                       + kxx * s_x_p * s_x_p //
                       + kxy * s_x_p * s_y_p //
                       + kxz * s_x_p * s_z_p //
                       + kxy * s_y_p * s_x_p //
                       + kyy * s_y_p * s_y_p //
                       + kyz * s_y_p * s_z_p //
                       + kxz * s_z_p * s_x_p //
                       + kyz * s_z_p * s_y_p //
                       + kzz * s_z_p * s_z_p;

      // Expansion

      const T div_s = (s_dsur_r + s_lsu) / s_sqrt_detq;

      const T kmm = 0                                //
                    + krr * (qurr - sur * sur)       //
                    + krt * (qurt - sur * sut)       //
                    + s_krp * (z_qurp - sur * z_sup) //
                    + krt * (qurt - sut * sur)       //
                    + ktt * (qutt - sut * sut)       //
                    + s_ktp * (z_qutp - sut * z_sup) //
                    + s_krp * (z_qurp - z_sup * sur) //
                    + s_ktp * (z_qutp - z_sup * sut) //
                    + ss_kpp * (zz_qupp - z_sup * z_sup);

      const T Theta = div_s - kmm;

      Thetaij(i, j) = Theta;
    }
  }

  const T area = sqrt(4 * T(M_PI)) * real(Alm(0, 0));

  const T cx = real(cxlm(0, 0)) / real(Alm(0, 0));
  const T cy = real(cylm(0, 0)) / real(Alm(0, 0));
  const T cz = real(czlm(0, 0)) / real(Alm(0, 0));

  const alm_t<T> Thetalm = expand(Thetaij);
  // [arXiv:gr-qc/0702038], (28)
  aij_t<T> Sij(geom);
  for (int i = 0; i < geom.ntheta; ++i)
#pragma omp simd
    for (int j = 0; j < geom.nphi; ++j)
      Sij(i, j) = lambdaij(i, j) * Thetaij(i, j);
  const alm_t<T> Slm = expand(Sij);

  alm_t<T> hlm_new(geom);
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
} // namespace AHFinder

////////////////////////////////////////////////////////////////////////////////

template <typename T>
expansion_t<T> update(const cGH *const cctkGH, const alm_t<T> &hlm) {
  const auto hij = evaluate(hlm);
  const auto coords = coords_from_shape(hij);
  // CCTK_VWARN(CCTK_WARN_ALERT, "Using Brill-Lindquist metric");
  // const auto metric = brill_lindquist_metric(cctkGH, coords);
  const auto metric = interpolate_metric(cctkGH, coords);
  const auto res = expansion(metric, hlm);
  return res;
}

template <typename T>
expansion_t<T> solve(const cGH *const cctkGH, const alm_t<T> &hlm_ini) {
  DECLARE_CCTK_PARAMETERS;
  int iter = 0;
  unique_ptr<const alm_t<T> > hlm_ptr =
      make_unique<alm_t<T> >(filter(hlm_ini, lmax_filter));
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

    alm_t<T> dhlm(geom);
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
    const auto dh2lm = expand(dh2ij);
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
    const auto Theta2lm = expand(Theta2ij);
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
    hlm_ptr = make_unique<alm_t<T> >(move(hlm_new));
  }
}

////////////////////////////////////////////////////////////////////////////////

extern "C" void AHFinder_find(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_find;
  DECLARE_CCTK_PARAMETERS;

  const geom_t geom(npoints);

  const auto hlm = coefficients_from_const(geom, r0, r1z);
  const auto res = solve(cctkGH, hlm);
}

} // namespace AHFinder
