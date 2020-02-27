// Coefficients have a "sin(theta)" weight: The power of sin(theta)
// that they contain. This sin(theta) power is not stored. For
// example, [xyz]_t have a weight of 0, while [xyz]_p have a weight of
// 1. The stored coefficients are thus never singular at the poles,
// i.e. neither infinite nor always zero.
//
// We use a prefix "s_" to indicate a weight of 1, and a prefix "z_"
// (an upside-down "s_") to indicate a weight of -1.
//
// Example: qpp = sin(theta)^2 * ss_qpp

#include <ssht.h>

#include <dual.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cassert>
#include <cmath>
#include <complex>
#include <iostream>
#include <vector>

namespace AHFinder {
using namespace std;

extern "C" void AHFinder_find(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_find;
  DECLARE_CCTK_PARAMETERS;

  // Rename parameters to avoid confusion
  const int filter_lmax = lmax;
  const int nmodes = npoints;
  {

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

    ////////////////////////////////////////////////////////////////////////////

    // Initial conditions for h^lm

    vector<CCTK_COMPLEX> hlm(ncoeffs);
    for (int l = 0; l <= lmax; ++l)
      for (int m = -l; m <= l; ++m)
        hlm.at(cind(l, m)) = l == 0 && m == 0 ? sqrt(4 * M_PI) * r0 : 0;

    // Filter h^lm
    for (int l = filter_lmax + 1; l <= lmax; ++l)
      for (int m = -l; m <= l; ++m)
        hlm.at(cind(l, m)) = 0;

    // Evaluate h^ij

    vector<CCTK_REAL> hij(npoints);
    ssht_core_mw_inverse_sov_sym_real(hij.data(), hlm.data(), nmodes, method,
                                      verbosity);

    // Evaluate derivatives of h^ij

    const int dh_spin = 1;
    vector<CCTK_COMPLEX> dhlm(ncoeffs);
    for (int l = 0; l <= lmax; ++l)
      for (int m = -l; m <= l; ++m)
        dhlm.at(cind(l, m)) = sqrt(CCTK_REAL(l * (l + 1))) * hlm.at(cind(l, m));

    vector<CCTK_COMPLEX> dhij(npoints);
    ssht_core_mw_inverse_sov_sym(dhij.data(), dhlm.data(), nmodes, dh_spin,
                                 method, verbosity);

    // Calculate interpolation coordinates

    vector<CCTK_REAL> coordsx(npoints), coordsy(npoints), coordsz(npoints);
    for (int i = 0; i < ntheta; ++i) {
      for (int j = 0; j < nphi; ++j) {
        const int ind2d = gind(i, j);
        const CCTK_REAL r = hij[ind2d];
        const CCTK_REAL theta = coord_theta(i, j);
        const CCTK_REAL phi = coord_phi(i, j);
        coordsx[ind2d] = r * sin(theta) * cos(phi);
        coordsy[ind2d] = r * sin(theta) * sin(phi);
        coordsz[ind2d] = r * cos(theta);
      }
    }

    // Interpolate metric and extrinsic curvature

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
    const vector<CCTK_INT> varinds{
        gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
        gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
        gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
        gxx_ind, gxy_ind, gxz_ind, gyy_ind, gyz_ind, gzz_ind, //
        kxx_ind, kxy_ind, kxz_ind, kyy_ind, kyz_ind, kzz_ind, //
    };
    const vector<CCTK_INT> operations{
        0, 0, 0, 0, 0, 0, //
        1, 1, 1, 1, 1, 1, //
        2, 2, 2, 2, 2, 2, //
        3, 3, 3, 3, 3, 3, //
        0, 0, 0, 0, 0, 0, //
    };
    constexpr int nvars = 6 * (1 + 3 + 1);
    assert(int(varinds.size()) == nvars);
    assert(int(operations.size()) == nvars);
    int index = 0;
    const int igxx = index++;
    const int igxy = index++;
    const int igxz = index++;
    const int igyy = index++;
    const int igyz = index++;
    const int igzz = index++;
    const int igxx_x = index++;
    const int igxy_x = index++;
    const int igxz_x = index++;
    const int igyy_x = index++;
    const int igyz_x = index++;
    const int igzz_x = index++;
    const int igxx_y = index++;
    const int igxy_y = index++;
    const int igxz_y = index++;
    const int igyy_y = index++;
    const int igyz_y = index++;
    const int igzz_y = index++;
    const int igxx_z = index++;
    const int igxy_z = index++;
    const int igxz_z = index++;
    const int igyy_z = index++;
    const int igyz_z = index++;
    const int igzz_z = index++;
    const int ikxx = index++;
    const int ikxy = index++;
    const int ikxz = index++;
    const int ikyy = index++;
    const int ikyz = index++;
    const int ikzz = index++;
    assert(index == nvars);

    vector<vector<CCTK_REAL> > gij(nvars);
    {
      vector<CCTK_REAL *> gijptrs(nvars);
      for (int n = 0; n < nvars; ++n) {
        gij.at(n).resize(npoints);
        gijptrs.at(n) = gij.at(n).data();
      }
      Interpolate(cctkGH, npoints, coordsx.data(), coordsy.data(),
                  coordsz.data(), nvars, varinds.data(), operations.data(),
                  gijptrs.data());
    }

    // Expand metric and extrinsic curvature

    vector<vector<CCTK_COMPLEX> > glm(nvars);
    for (int n = 0; n < nvars; ++n) {
      glm.at(n).resize(ncoeffs);
      ssht_core_mw_forward_sov_conv_sym_real(glm.at(n).data(), gij.at(n).data(),
                                             nmodes, method, verbosity);
    }

    // Filter metric and extrinsic curvature

    for (int n = 0; n < nvars; ++n)
      for (int l = filter_lmax + 1; l <= lmax; ++l)
        for (int m = -l; m <= l; ++m)
          glm.at(n).at(cind(l, m)) = 0;

    // Evaluate metric and extrinsic curvature

    for (int n = 0; n < nvars; ++n)
      ssht_core_mw_inverse_sov_sym_real(gij.at(n).data(), glm.at(n).data(),
                                        nmodes, method, verbosity);

    // Calculate spacelike normal s^i and its radial derivative s^i,r

    vector<CCTK_REAL> sxij(npoints), syij(npoints), szij(npoints);
    vector<CCTK_REAL> sqrtdetg_sxij(npoints), sqrtdetg_syij(npoints),
        sqrtdetg_szij(npoints);
    vector<CCTK_REAL> sqrtdetg_sx_rij(npoints), sqrtdetg_sy_rij(npoints),
        sqrtdetg_sz_rij(npoints);

    for (int i = 0; i < ntheta; ++i) {
      for (int j = 0; j < nphi; ++j) {
        const int ind2d = gind(i, j);

        // Coordinates

        // Cartesian: x y z
        // spherical: R Theta Phi
        // horizon:   r theta phi

        // x = R * sin(Theta) * cos(Phi)
        // y = R * sin(Theta) * sin(Phi)
        // z = R * cos(Theta)

        // r     = R - h(Theta, Phi)
        // theta = Theta
        // phi   = Phi

        const CCTK_REAL h = hij[ind2d];
        // dX = - \dh X = - (d/dTheta X + i/sin Theta d/dPhi X)
        const CCTK_COMPLEX dh = dhij[ind2d];
        const CCTK_REAL h_T = -real(dh);
        const CCTK_REAL s_h_P = -imag(dh);

        const CCTK_REAL R = h;
        const CCTK_REAL Theta = coord_theta(i, j);
        const CCTK_REAL Phi = coord_phi(i, j);

        // const CCTK_REAL r = R - h;
        // const CCTK_REAL theta = Theta;
        // const CCTK_REAL phi = Phi;

        // only non-zero terms
        const CCTK_REAL r_R = 1;
        const CCTK_REAL r_T = -h_T;
        const CCTK_REAL s_r_P = -s_h_P;
        const CCTK_REAL t_T = 1;
        const CCTK_REAL p_P = 1;

        // inverse
        const CCTK_REAL R_r = r_R;
        const CCTK_REAL R_t = r_T;
        const CCTK_REAL s_R_p = s_r_P;
        const CCTK_REAL T_t = t_T;
        const CCTK_REAL P_p = p_P;

        // x^i
        // const CCTK_REAL x = R * sin(Theta) * cos(Phi);
        // const CCTK_REAL y = R * sin(Theta) * sin(Phi);
        // const CCTK_REAL z = R * cos(Theta);

        // only non-zero terms
        const CCTK_REAL x_R = sin(Theta) * cos(Phi);
        const CCTK_REAL x_T = R * cos(Theta) * cos(Phi);
        const CCTK_REAL s_x_P = -R * sin(Phi);
        const CCTK_REAL y_R = sin(Theta) * sin(Phi);
        const CCTK_REAL y_T = R * cos(Theta) * sin(Phi);
        const CCTK_REAL s_y_P = R * cos(Phi);
        const CCTK_REAL z_R = cos(Theta);
        const CCTK_REAL z_T = -R * sin(Theta);

        // only radial derivatives
        const CCTK_REAL x_RT = cos(Theta) * cos(Phi);
        const CCTK_REAL s_x_RP = -sin(Phi);
        const CCTK_REAL y_RT = cos(Theta) * sin(Phi);
        const CCTK_REAL s_y_RP = cos(Phi);
        const CCTK_REAL z_RT = -sin(Theta);

        const CCTK_REAL x_r = x_R * R_r;
        const CCTK_REAL x_t = x_R * R_t + x_T * T_t;
        const CCTK_REAL s_x_p = x_R * s_R_p + s_x_P * P_p;
        const CCTK_REAL y_r = y_R * R_r;
        const CCTK_REAL y_t = y_R * R_t + y_T * T_t;
        const CCTK_REAL s_y_p = y_R * s_R_p + s_y_P * P_p;
        const CCTK_REAL z_r = z_R * R_r;
        const CCTK_REAL z_t = z_R * R_t + z_T * T_t;
        const CCTK_REAL s_z_p = z_R * s_R_p;

        // only radial derivatives
        const CCTK_REAL x_rr = 0;
        const CCTK_REAL x_rt = x_RT * R_r * T_t;
        const CCTK_REAL s_x_rp = s_x_RP * R_r * P_p;
        const CCTK_REAL y_rr = 0;
        const CCTK_REAL y_rt = y_RT * R_r * T_t;
        const CCTK_REAL s_y_rp = s_y_RP * R_r * P_p;
        const CCTK_REAL z_rr = 0;
        const CCTK_REAL z_rt = z_RT * R_r * T_t;
        const CCTK_REAL s_z_rp = 0;

        // Three-metric
        const CCTK_REAL gxx0 = gij.at(igxx)[ind2d];
        const CCTK_REAL gxy0 = gij.at(igxy)[ind2d];
        const CCTK_REAL gxz0 = gij.at(igxz)[ind2d];
        const CCTK_REAL gyy0 = gij.at(igyy)[ind2d];
        const CCTK_REAL gyz0 = gij.at(igyz)[ind2d];
        const CCTK_REAL gzz0 = gij.at(igzz)[ind2d];

        const CCTK_REAL gxx_x = gij.at(igxx_x)[ind2d];
        const CCTK_REAL gxy_x = gij.at(igxy_x)[ind2d];
        const CCTK_REAL gxz_x = gij.at(igxz_x)[ind2d];
        const CCTK_REAL gyy_x = gij.at(igyy_x)[ind2d];
        const CCTK_REAL gyz_x = gij.at(igyz_x)[ind2d];
        const CCTK_REAL gzz_x = gij.at(igzz_x)[ind2d];
        const CCTK_REAL gxx_y = gij.at(igxx_y)[ind2d];
        const CCTK_REAL gxy_y = gij.at(igxy_y)[ind2d];
        const CCTK_REAL gxz_y = gij.at(igxz_y)[ind2d];
        const CCTK_REAL gyy_y = gij.at(igyy_y)[ind2d];
        const CCTK_REAL gyz_y = gij.at(igyz_y)[ind2d];
        const CCTK_REAL gzz_y = gij.at(igzz_y)[ind2d];
        const CCTK_REAL gxx_z = gij.at(igxx_z)[ind2d];
        const CCTK_REAL gxy_z = gij.at(igxy_z)[ind2d];
        const CCTK_REAL gxz_z = gij.at(igxz_z)[ind2d];
        const CCTK_REAL gyy_z = gij.at(igyy_z)[ind2d];
        const CCTK_REAL gyz_z = gij.at(igyz_z)[ind2d];
        const CCTK_REAL gzz_z = gij.at(igzz_z)[ind2d];

        // Radial derivative of metric
        const CCTK_REAL gxx_r = gxx_x * x_r + gxx_y * y_r + gxx_z * z_r;
        const CCTK_REAL gxy_r = gxy_x * x_r + gxy_y * y_r + gxy_z * z_r;
        const CCTK_REAL gxz_r = gxz_x * x_r + gxz_y * y_r + gxz_z * z_r;
        const CCTK_REAL gyy_r = gyy_x * x_r + gyy_y * y_r + gyy_z * z_r;
        const CCTK_REAL gyz_r = gyz_x * x_r + gyz_y * y_r + gyz_z * z_r;
        const CCTK_REAL gzz_r = gzz_x * x_r + gzz_y * y_r + gzz_z * z_r;

        // Dual quantities to handle radial derivatives

        const dual<CCTK_REAL> gxx{gxx0, gxx_r};
        const dual<CCTK_REAL> gxy{gxy0, gxy_r};
        const dual<CCTK_REAL> gxz{gxz0, gxz_r};
        const dual<CCTK_REAL> gyy{gyy0, gyy_r};
        const dual<CCTK_REAL> gyz{gyz0, gyz_r};
        const dual<CCTK_REAL> gzz{gzz0, gzz_r};

        // Determinant of three-metric
        const dual<CCTK_REAL> detg = -pow2(gxz) * gyy + 2 * gxy * gxz * gyz -
                                     gxx * pow2(gyz) - pow2(gxy) * gzz +
                                     gxx * gyy * gzz;

        // Triad m1 m2 s
        // m1^2 = m2^2 = s^2 = 1
        // m1 m2 = m1 s = m2 s = 0
        dual<CCTK_REAL> m1x{x_t, x_rt};
        dual<CCTK_REAL> m1y{y_t, y_rt};
        dual<CCTK_REAL> m1z{z_t, z_rt};
        dual<CCTK_REAL> m1lx = gxx * m1x + gxy * m1y + gxz * m1z;
        dual<CCTK_REAL> m1ly = gxy * m1x + gyy * m1y + gyz * m1z;
        dual<CCTK_REAL> m1lz = gxz * m1x + gyz * m1y + gzz * m1z;
        dual<CCTK_REAL> m12 = m1lx * m1x + m1ly * m1y + m1lz * m1z;
        m1x /= sqrt(m12);
        m1y /= sqrt(m12);
        m1z /= sqrt(m12);
        m1lx = gxx * m1x + gxy * m1y + gxz * m1z;
        m1ly = gxy * m1x + gyy * m1y + gyz * m1z;
        m1lz = gxz * m1x + gyz * m1y + gzz * m1z;
#ifdef CCTK_DEBUG
        // Test
        m12 = m1lx * m1x + m1ly * m1y + m1lz * m1z;
        assert(fabs(m12 - 1) <= 1.0e-12);
#endif

        dual<CCTK_REAL> m2x{s_x_p, s_x_rp};
        dual<CCTK_REAL> m2y{s_y_p, s_y_rp};
        dual<CCTK_REAL> m2z{s_z_p, s_z_rp};
        dual<CCTK_REAL> m1m2 = m1lx * m2x + m1ly * m2y + m1lz * m2z;
        m2x -= m1m2 * m1x;
        m2y -= m1m2 * m1y;
        m2z -= m1m2 * m1z;
        dual<CCTK_REAL> m2lx = gxx * m2x + gxy * m2y + gxz * m2z;
        dual<CCTK_REAL> m2ly = gxy * m2x + gyy * m2y + gyz * m2z;
        dual<CCTK_REAL> m2lz = gxz * m2x + gyz * m2y + gzz * m2z;
        dual<CCTK_REAL> m22 = m2lx * m2x + m2ly * m2y + m2lz * m2z;
        m2x /= sqrt(m22);
        m2y /= sqrt(m22);
        m2z /= sqrt(m22);
        m2lx = gxx * m2x + gxy * m2y + gxz * m2z;
        m2ly = gxy * m2x + gyy * m2y + gyz * m2z;
        m2lz = gxz * m2x + gyz * m2y + gzz * m2z;
#ifdef CCTK_DEBUG
        // Test
        m1m2 = m1lx * m2x + m1ly * m2y + m1lz * m2z;
        assert(fabs(m1m2) <= 1.0e-12);
        m22 = m2lx * m2x + m2ly * m2y + m2lz * m2z;
        assert(fabs(m22 - 1) <= 1.0e-12);
#endif

        dual<CCTK_REAL> sx{x_r, x_rr};
        dual<CCTK_REAL> sy{y_r, y_rr};
        dual<CCTK_REAL> sz{z_r, z_rr};
        dual<CCTK_REAL> m1s = m1lx * sx + m1ly * sy + m1lz * sz;
        sx -= m1s * m1x;
        sy -= m1s * m1y;
        sz -= m1s * m1z;
        dual<CCTK_REAL> m2s = m2lx * sx + m2ly * sy + m2lz * sz;
        sx -= m2s * m2x;
        sy -= m2s * m2y;
        sz -= m2s * m2z;
        dual<CCTK_REAL> slx = gxx * sx + gxy * sy + gxz * sz;
        dual<CCTK_REAL> sly = gxy * sx + gyy * sy + gyz * sz;
        dual<CCTK_REAL> slz = gxz * sx + gyz * sy + gzz * sz;
        dual<CCTK_REAL> s2 = slx * sx + sly * sy + slz * sz;
        sx /= sqrt(s2);
        sy /= sqrt(s2);
        sz /= sqrt(s2);
        slx = gxx * sx + gxy * sy + gxz * sz;
        sly = gxy * sx + gyy * sy + gyz * sz;
        slz = gxz * sx + gyz * sy + gzz * sz;
#ifdef CCTK_DEBUG
        // Test
        m1s = m1lx * sx + m1ly * sy + m1lz * sz;
        assert(fabs(m1s) <= 1.0e-12);
        m2s = m2lx * sx + m2ly * sy + m2lz * sz;
        assert(fabs(m2s) <= 1.0e-12);
        s2 = slx * sx + sly * sy + slz * sz;
        assert(fabs(s2 - 1) <= 1.0e-12);
#endif

        // Store spacelike normal
        sxij[ind2d] = sx.val;
        syij[ind2d] = sy.val;
        szij[ind2d] = sz.val;

        sqrtdetg_sxij[ind2d] = (sqrt(detg) * sx).val;
        sqrtdetg_syij[ind2d] = (sqrt(detg) * sy).val;
        sqrtdetg_szij[ind2d] = (sqrt(detg) * sz).val;

        sqrtdetg_sx_rij[ind2d] = (sqrt(detg) * sx).eps;
        sqrtdetg_sy_rij[ind2d] = (sqrt(detg) * sy).eps;
        sqrtdetg_sz_rij[ind2d] = (sqrt(detg) * sz).eps;
      }
    }

    // Calculate divergence of spacelike normal

    vector<CCTK_COMPLEX> sqrtdetg_sxlm(npoints), sqrtdetg_sylm(npoints),
        sqrtdetg_szlm(npoints);
    ssht_core_mw_forward_sov_conv_sym_real(
        sqrtdetg_sxlm.data(), sqrtdetg_sxij.data(), nmodes, method, verbosity);
    ssht_core_mw_forward_sov_conv_sym_real(
        sqrtdetg_sylm.data(), sqrtdetg_syij.data(), nmodes, method, verbosity);
    ssht_core_mw_forward_sov_conv_sym_real(
        sqrtdetg_szlm.data(), sqrtdetg_szij.data(), nmodes, method, verbosity);

    vector<CCTK_COMPLEX> d_sqrtdetg_sxlm(npoints), d_sqrtdetg_sylm(npoints),
        d_sqrtdetg_szlm(npoints);
    for (int l = 0; l <= lmax; ++l) {
      for (int m = -l; m <= l; ++m) {
        d_sqrtdetg_sxlm.at(cind(l, m)) =
            sqrt(CCTK_REAL(l * (l + 1))) * sqrtdetg_sxlm.at(cind(l, m));
        d_sqrtdetg_sylm.at(cind(l, m)) =
            sqrt(CCTK_REAL(l * (l + 1))) * sqrtdetg_sylm.at(cind(l, m));
        d_sqrtdetg_szlm.at(cind(l, m)) =
            sqrt(CCTK_REAL(l * (l + 1))) * sqrtdetg_szlm.at(cind(l, m));
      }
    }

    const int ds_spin = 1;
    vector<CCTK_COMPLEX> d_sqrtdetg_sxij(npoints), d_sqrtdetg_syij(npoints),
        d_sqrtdetg_szij(npoints);
    ssht_core_mw_inverse_sov_sym(d_sqrtdetg_sxij.data(), d_sqrtdetg_sxlm.data(),
                                 nmodes, ds_spin, method, verbosity);
    ssht_core_mw_inverse_sov_sym(d_sqrtdetg_syij.data(), d_sqrtdetg_sylm.data(),
                                 nmodes, ds_spin, method, verbosity);
    ssht_core_mw_inverse_sov_sym(d_sqrtdetg_szij.data(), d_sqrtdetg_szlm.data(),
                                 nmodes, ds_spin, method, verbosity);

    // Calculate expansion

    vector<CCTK_REAL> Theta_ij(npoints);

    for (int i = 0; i < ntheta; ++i) {
      for (int j = 0; j < nphi; ++j) {
        const int ind2d = gind(i, j);

        // Coordinates

        const CCTK_REAL h = hij[ind2d];
        // dX = \dh X = - (d/dTheta X + i/sin Theta d/dPhi X)
        const CCTK_COMPLEX dh = dhij[ind2d];
        const CCTK_REAL h_T = -real(dh);
        const CCTK_REAL s_h_P = -imag(dh);

        const CCTK_REAL R = h;
        const CCTK_REAL Theta = coord_theta(i, j);
        const CCTK_REAL Phi = coord_phi(i, j);

        // const CCTK_REAL r = R - h;
        // const CCTK_REAL theta = Theta;
        // const CCTK_REAL phi = Phi;

        // only non-zero terms
        const CCTK_REAL r_R = 1;
        const CCTK_REAL r_T = -h_T;
        const CCTK_REAL s_r_P = -s_h_P;
        const CCTK_REAL t_T = 1;
        const CCTK_REAL p_P = 1;

        // inverse
        const CCTK_REAL R_r = r_R;
        const CCTK_REAL R_t = r_T;
        const CCTK_REAL s_R_p = s_r_P;
        const CCTK_REAL T_t = t_T;
        const CCTK_REAL P_p = p_P;

        // x^i
        // const CCTK_REAL x = R * sin(Theta) * cos(Phi);
        // const CCTK_REAL y = R * sin(Theta) * sin(Phi);
        // const CCTK_REAL z = R * cos(Theta);

        // only non-zero terms
        const CCTK_REAL x_R = sin(Theta) * cos(Phi);
        const CCTK_REAL x_T = R * cos(Theta) * cos(Phi);
        const CCTK_REAL s_x_P = -R * sin(Phi);
        const CCTK_REAL y_R = sin(Theta) * sin(Phi);
        const CCTK_REAL y_T = R * cos(Theta) * sin(Phi);
        const CCTK_REAL s_y_P = R * cos(Phi);
        const CCTK_REAL z_R = cos(Theta);
        const CCTK_REAL z_T = -R * sin(Theta);

        const CCTK_REAL x_r = x_R * R_r;
        const CCTK_REAL x_t = x_R * R_t + x_T * T_t;
        const CCTK_REAL s_x_p = x_R * s_R_p + s_x_P * P_p;
        const CCTK_REAL y_r = y_R * R_r;
        const CCTK_REAL y_t = y_R * R_t + y_T * T_t;
        const CCTK_REAL s_y_p = y_R * s_R_p + s_y_P * P_p;
        const CCTK_REAL z_r = z_R * R_r;
        const CCTK_REAL z_t = z_R * R_t + z_T * T_t;
        const CCTK_REAL s_z_p = z_R * s_R_p;

        const CCTK_REAL s_det_xr = -s_z_p * x_t * y_r + s_z_p * x_r * y_t +
                                   s_y_p * x_t * z_r - s_x_p * y_t * z_r -
                                   s_y_p * x_r * z_t + s_x_p * y_r * z_t;

        const CCTK_REAL r_x = (s_z_p * y_t - s_y_p * z_t) / s_det_xr;
        const CCTK_REAL r_y = (-s_z_p * x_t + s_x_p * z_t) / s_det_xr;
        const CCTK_REAL r_z = (s_y_p * x_t - s_x_p * y_t) / s_det_xr;
        const CCTK_REAL theta_x = (-s_z_p * y_r + s_y_p * z_r) / s_det_xr;
        const CCTK_REAL theta_y = (s_z_p * x_r - s_x_p * z_r) / s_det_xr;
        const CCTK_REAL theta_z = (-s_y_p * x_r + s_x_p * y_r) / s_det_xr;
        const CCTK_REAL z_phi_x = (-y_t * z_r + y_r * z_t) / s_det_xr;
        const CCTK_REAL z_phi_y = (x_t * z_r - x_r * z_t) / s_det_xr;
        const CCTK_REAL z_phi_z = (-x_t * y_r + x_r * y_t) / s_det_xr;
#ifdef CCTK_DEBUG
        // Test
        assert(fabs(r_x * x_r + r_y * y_r + r_z * z_r - 1) <= 1.0e-12);
        assert(fabs(r_x * x_t + r_y * y_t + r_z * z_t) <= 1.0e-12);
        assert(fabs(r_x * s_x_p + r_y * s_y_p + r_z * s_z_p) <= 1.0e-12);
        assert(fabs(theta_x * x_r + theta_y * y_r + theta_z * z_r) <= 1.0e-12);
        assert(fabs(theta_x * x_t + theta_y * y_t + theta_z * z_t - 1) <=
               1.0e-12);
        assert(fabs(theta_x * s_x_p + theta_y * s_y_p + theta_z * s_z_p) <=
               1.0e-12);
        assert(fabs(z_phi_x * x_r + z_phi_y * y_r + z_phi_z * z_r) <= 1.0e-12);
        assert(fabs(z_phi_x * x_t + z_phi_y * y_t + z_phi_z * z_t) <= 1.0e-12);
        assert(fabs(z_phi_x * s_x_p + z_phi_y * s_y_p + z_phi_z * s_z_p - 1) <=
               1.0e-12);
#endif

        // Read three-metric and extrinsic curvature

        const CCTK_REAL gxx = gij.at(igxx)[ind2d];
        const CCTK_REAL gxy = gij.at(igxy)[ind2d];
        const CCTK_REAL gxz = gij.at(igxz)[ind2d];
        const CCTK_REAL gyy = gij.at(igyy)[ind2d];
        const CCTK_REAL gyz = gij.at(igyz)[ind2d];
        const CCTK_REAL gzz = gij.at(igzz)[ind2d];
        const CCTK_REAL kxx = gij.at(ikxx)[ind2d];
        const CCTK_REAL kxy = gij.at(ikxy)[ind2d];
        const CCTK_REAL kxz = gij.at(ikxz)[ind2d];
        const CCTK_REAL kyy = gij.at(ikyy)[ind2d];
        const CCTK_REAL kyz = gij.at(ikyz)[ind2d];
        const CCTK_REAL kzz = gij.at(ikzz)[ind2d];

        // Determinant of three-metric
        const CCTK_REAL detg = -pow(gxz, 2) * gyy + 2 * gxy * gxz * gyz -
                               gxx * pow(gyz, 2) - pow(gxy, 2) * gzz +
                               gxx * gyy * gzz;

        // Inverse three-metric
        const CCTK_REAL guxx = (-pow(gyz, 2) + gyy * gzz) / detg;
        const CCTK_REAL guxy = (gxz * gyz - gxy * gzz) / detg;
        const CCTK_REAL guxz = (-gxz * gyy + gxy * gyz) / detg;
        const CCTK_REAL guyy = (-pow(gxz, 2) + gxx * gzz) / detg;
        const CCTK_REAL guyz = (gxy * gxz - gxx * gyz) / detg;
        const CCTK_REAL guzz = (-pow(gxy, 2) + gxx * gyy) / detg;

        // Read spacelike normal

        const CCTK_REAL sx = sxij[ind2d];
        const CCTK_REAL sy = syij[ind2d];
        const CCTK_REAL sz = szij[ind2d];

        // Read derivatives of densitized spacelike normal

        const CCTK_REAL sqrtdetg_sx_r = sqrtdetg_sx_rij[ind2d];
        const CCTK_REAL sqrtdetg_sy_r = sqrtdetg_sy_rij[ind2d];
        const CCTK_REAL sqrtdetg_sz_r = sqrtdetg_sz_rij[ind2d];

        const CCTK_COMPLEX d_sqrtdetg_sx = d_sqrtdetg_sxij[ind2d];
        const CCTK_COMPLEX d_sqrtdetg_sy = d_sqrtdetg_syij[ind2d];
        const CCTK_COMPLEX d_sqrtdetg_sz = d_sqrtdetg_szij[ind2d];

        const CCTK_REAL sqrtdetg_sx_t = -real(d_sqrtdetg_sx);
        const CCTK_REAL sqrtdetg_sy_t = -real(d_sqrtdetg_sy);
        const CCTK_REAL sqrtdetg_sz_t = -real(d_sqrtdetg_sz);

        const CCTK_REAL s_sqrtdetg_sx_p = -imag(d_sqrtdetg_sx);
        const CCTK_REAL s_sqrtdetg_sy_p = -imag(d_sqrtdetg_sy);
        const CCTK_REAL s_sqrtdetg_sz_p = -imag(d_sqrtdetg_sz);

        const CCTK_REAL sqrtdetg_sx_x = sqrtdetg_sx_r * r_x +
                                        sqrtdetg_sx_t * theta_x +
                                        s_sqrtdetg_sx_p * z_phi_x;
        const CCTK_REAL sqrtdetg_sy_y = sqrtdetg_sy_r * r_y +
                                        sqrtdetg_sy_t * theta_y +
                                        s_sqrtdetg_sy_p * z_phi_y;
        const CCTK_REAL sqrtdetg_sz_z = sqrtdetg_sz_r * r_z +
                                        sqrtdetg_sz_t * theta_z +
                                        s_sqrtdetg_sz_p * z_phi_z;
        const CCTK_REAL div_sqrtdetg_s =
            sqrtdetg_sx_x + sqrtdetg_sy_y + sqrtdetg_sz_z;
        const CCTK_REAL div_s = div_sqrtdetg_s / sqrt(detg);

        // Theta_(l) = D_i s^i + K_ij s^i s^j - K   [arxiv:gr-qc/0512169] (3.1)
        const CCTK_REAL Theta_l = div_s                    //
                                  + kxx * (sx * sx - guxx) //
                                  + kxy * (sx * sy - guxy) //
                                  + kxz * (sx * sz - guxz) //
                                  + kxy * (sy * sx - guxy) //
                                  + kyy * (sy * sy - guyy) //
                                  + kyz * (sy * sz - guyz) //
                                  + kxz * (sz * sx - guxz) //
                                  + kyz * (sz * sy - guyz) //
                                  + kzz * (sz * sz - guzz);

        // Store expansion
        Theta_ij[ind2d] = Theta_l;
      }
    }

    double h_min = INFINITY, h_max = -INFINITY;
    double Theta_min = INFINITY, Theta_max = -INFINITY;
    for (int i = 0; i < ntheta; ++i) {
      for (int j = 0; j < nphi; ++j) {
        const int ind2d = gind(i, j);
        const CCTK_REAL h = hij[ind2d];
        const CCTK_REAL Theta = Theta_ij[ind2d];
        h_min = fmin(h_min, h);
        h_max = fmax(h_max, h);
        Theta_min = fmin(Theta_min, Theta);
        Theta_max = fmax(Theta_max, Theta);
      }
    }
    cout << "h in [" << h_min << "; " << h_max << "]\n";
    cout << "Theta in [" << Theta_min << "; " << Theta_max << "]\n";
  }
} // namespace AHFinder

} // namespace AHFinder
