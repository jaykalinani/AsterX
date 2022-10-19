#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include "utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

extern "C" void AsterX_SourceTerms(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SourceTerms;
  DECLARE_CCTK_PARAMETERS;

  if ((local_spatial_order != 2) && (local_spatial_order != 4)) {
    CCTK_VERROR("local_spatial_order must be set to 2 or 4.");
  }

  /* Loop over the entire grid (0 to n-1 cells in each direction) */
  grid.loop_all_device<
      1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    /* TODO: local_spatial_order = spatial_order; //(taken from ADMMacros) */
    /* order of finite differencing */
    const CCTK_INT local_spatial_order = 2;

    /* Computing metric components at cell centers */
    CCTK_REAL alp_avg = calc_avg_v2c(alp, p);
    CCTK_REAL betax_avg = calc_avg_v2c(betax, p);
    CCTK_REAL betay_avg = calc_avg_v2c(betay, p);
    CCTK_REAL betaz_avg = calc_avg_v2c(betaz, p);

    CCTK_REAL gxx_avg = calc_avg_v2c(gxx, p);
    CCTK_REAL gxy_avg = calc_avg_v2c(gxy, p);
    CCTK_REAL gxz_avg = calc_avg_v2c(gxz, p);
    CCTK_REAL gyy_avg = calc_avg_v2c(gyy, p);
    CCTK_REAL gyz_avg = calc_avg_v2c(gyz, p);
    CCTK_REAL gzz_avg = calc_avg_v2c(gzz, p);

    CCTK_REAL kxx_avg = calc_avg_v2c(kxx, p);
    CCTK_REAL kxy_avg = calc_avg_v2c(kxy, p);
    CCTK_REAL kxz_avg = calc_avg_v2c(kxz, p);
    CCTK_REAL kyy_avg = calc_avg_v2c(kyy, p);
    CCTK_REAL kyz_avg = calc_avg_v2c(kyz, p);
    CCTK_REAL kzz_avg = calc_avg_v2c(kzz, p);

    /* Determinant of spatial metric */
    const smat<CCTK_REAL, 3, DN, DN> g{gxx_avg, gxy_avg, gxz_avg,
                                       gyy_avg, gyz_avg, gzz_avg};
    const CCTK_REAL detg = calc_det(g);

    const CCTK_REAL sqrt_detg = sqrt(detg);

    /* Upper metric */
    const array<CCTK_REAL, 6> ug =
        calc_upperg(gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg, detg);

    /* Computing v_j */
    const array<CCTK_REAL, 3> v_up = {velx(p.I), vely(p.I), velz(p.I)};
    const array<CCTK_REAL, 3> v_low =
        calc_vlow(v_up, gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg);

    const CCTK_REAL w_lorentz = calc_wlor(v_low, v_up);
    /* Computing [ \rho(1+\epsilon) + Pgas ]*W^2 */
    const CCTK_REAL rhoenthalpyW2 =
        (rho(p.I) * (1.0 + eps(p.I)) + press(p.I)) * w_lorentz * w_lorentz;

    /* TODO: define shift_state: storage for shift alloted if shift_state = 1 */
    /* For now, we assume shift has storage */

    // if (shift_state != 0){

    //    CCTK_REAL shiftx = betax(i,j,k)
    //    shifty = betay(i,j,k)
    //    shiftz = betaz(i,j,k)

    /* Derivatives of the lapse, shift and metric */

    CCTK_REAL dx_alp, dy_alp, dz_alp;
    CCTK_REAL dx_betax, dy_betax, dz_betax;
    CCTK_REAL dx_betay, dy_betay, dz_betay;
    CCTK_REAL dx_betaz, dy_betaz, dz_betaz;
    CCTK_REAL dx_gxx, dy_gxx, dz_gxx;
    CCTK_REAL dx_gxy, dy_gxy, dz_gxy;
    CCTK_REAL dx_gxz, dy_gxz, dz_gxz;
    CCTK_REAL dx_gyy, dy_gyy, dz_gyy;
    CCTK_REAL dx_gyz, dy_gyz, dz_gyz;
    CCTK_REAL dx_gzz, dy_gzz, dz_gzz;

    if (local_spatial_order == 2) {
      /* calc_fd2_v2c takes vertex center input, computes edge-center
       * derivatives along either dir=0, 1 or 2 using 2nd finite difference,
       * i.e, at the four edges of the cube with p.I as starting point. The four
       * edge-centered values are then interpolated to the cell-center
       * using 2nd order interpolation */

      dx_alp = calc_fd2_v2c(alp, p, 0);
      dy_alp = calc_fd2_v2c(alp, p, 1);
      dz_alp = calc_fd2_v2c(alp, p, 2);

      dx_betax = calc_fd2_v2c(betax, p, 0);
      dy_betax = calc_fd2_v2c(betax, p, 1);
      dz_betax = calc_fd2_v2c(betax, p, 2);

      dx_betay = calc_fd2_v2c(betay, p, 0);
      dy_betay = calc_fd2_v2c(betay, p, 1);
      dz_betay = calc_fd2_v2c(betay, p, 2);

      dx_betaz = calc_fd2_v2c(betaz, p, 0);
      dy_betaz = calc_fd2_v2c(betaz, p, 1);
      dz_betaz = calc_fd2_v2c(betaz, p, 2);

      dx_gxx = calc_fd2_v2c(gxx, p, 0);
      dy_gxx = calc_fd2_v2c(gxx, p, 1);
      dz_gxx = calc_fd2_v2c(gxx, p, 2);

      dx_gxy = calc_fd2_v2c(gxy, p, 0);
      dy_gxy = calc_fd2_v2c(gxy, p, 1);
      dz_gxy = calc_fd2_v2c(gxy, p, 2);

      dx_gxz = calc_fd2_v2c(gxz, p, 0);
      dy_gxz = calc_fd2_v2c(gxz, p, 1);
      dz_gxz = calc_fd2_v2c(gxz, p, 2);

      dx_gyy = calc_fd2_v2c(gyy, p, 0);
      dy_gyy = calc_fd2_v2c(gyy, p, 1);
      dz_gyy = calc_fd2_v2c(gyy, p, 2);

      dx_gyz = calc_fd2_v2c(gyz, p, 0);
      dy_gyz = calc_fd2_v2c(gyz, p, 1);
      dz_gyz = calc_fd2_v2c(gyz, p, 2);

      dx_gzz = calc_fd2_v2c(gzz, p, 0);
      dy_gzz = calc_fd2_v2c(gzz, p, 1);
      dz_gzz = calc_fd2_v2c(gzz, p, 2);

    } else if (local_spatial_order == 4) {
      /* calc_fd4_v2c takes vertex center input, computes edge-center
       * derivatives along either dir=0, 1 or 2 using 4th order finite
       * difference/ The eight edge-centered values are then interpolated to the
       * center using 4th order interpolation */

      dx_alp = calc_fd4_v2c(alp, p, 0);
      dy_alp = calc_fd4_v2c(alp, p, 1);
      dz_alp = calc_fd4_v2c(alp, p, 2);

      dx_betax = calc_fd4_v2c(betax, p, 0);
      dy_betax = calc_fd4_v2c(betax, p, 1);
      dz_betax = calc_fd4_v2c(betax, p, 2);

      dx_betay = calc_fd4_v2c(betay, p, 0);
      dy_betay = calc_fd4_v2c(betay, p, 1);
      dz_betay = calc_fd4_v2c(betay, p, 2);

      dx_betaz = calc_fd4_v2c(betaz, p, 0);
      dy_betaz = calc_fd4_v2c(betaz, p, 1);
      dz_betaz = calc_fd4_v2c(betaz, p, 2);

      dx_gxx = calc_fd4_v2c(gxx, p, 0);
      dy_gxx = calc_fd4_v2c(gxx, p, 1);
      dz_gxx = calc_fd4_v2c(gxx, p, 2);

      dx_gxy = calc_fd4_v2c(gxy, p, 0);
      dy_gxy = calc_fd4_v2c(gxy, p, 1);
      dz_gxy = calc_fd4_v2c(gxy, p, 2);

      dx_gxz = calc_fd4_v2c(gxz, p, 0);
      dy_gxz = calc_fd4_v2c(gxz, p, 1);
      dz_gxz = calc_fd4_v2c(gxz, p, 2);

      dx_gyy = calc_fd4_v2c(gyy, p, 0);
      dy_gyy = calc_fd4_v2c(gyy, p, 1);
      dz_gyy = calc_fd4_v2c(gyy, p, 2);

      dx_gyz = calc_fd4_v2c(gyz, p, 0);
      dy_gyz = calc_fd4_v2c(gyz, p, 1);
      dz_gyz = calc_fd4_v2c(gyz, p, 2);

      dx_gzz = calc_fd4_v2c(gzz, p, 0);
      dy_gzz = calc_fd4_v2c(gzz, p, 1);
      dz_gzz = calc_fd4_v2c(gzz, p, 2);
    }

    const CCTK_REAL velxshift = v_up[0] - betax_avg / alp_avg;
    const CCTK_REAL velyshift = v_up[1] - betay_avg / alp_avg;
    const CCTK_REAL velzshift = v_up[2] - betaz_avg / alp_avg;

    /* Computing T_munu */
    CCTK_REAL t00 = (rhoenthalpyW2 - press(p.I)) / pow2(alp_avg);
    CCTK_REAL t0x = rhoenthalpyW2 * velxshift / alp_avg +
                    press(p.I) * betax_avg / pow2(alp_avg);
    CCTK_REAL t0y = rhoenthalpyW2 * velyshift / alp_avg +
                    press(p.I) * betay_avg / pow2(alp_avg);
    CCTK_REAL t0z = rhoenthalpyW2 * velzshift / alp_avg +
                    press(p.I) * betaz_avg / pow2(alp_avg);
    CCTK_REAL txx = rhoenthalpyW2 * velxshift * velxshift +
                    press(p.I) * (ug[0] - betax_avg * betax_avg /
                                              pow2(alp_avg)); // ug[0]=uxx
    CCTK_REAL txy = rhoenthalpyW2 * velxshift * velyshift +
                    press(p.I) * (ug[1] - betax_avg * betay_avg /
                                              pow2(alp_avg)); // ug[1]=uxy
    CCTK_REAL txz = rhoenthalpyW2 * velxshift * velzshift +
                    press(p.I) * (ug[2] - betax_avg * betaz_avg /
                                              pow2(alp_avg)); // ug[2]=uxz
    CCTK_REAL tyy = rhoenthalpyW2 * velyshift * velyshift +
                    press(p.I) * (ug[3] - betay_avg * betay_avg /
                                              pow2(alp_avg)); // ug[3]=uyy
    CCTK_REAL tyz = rhoenthalpyW2 * velyshift * velzshift +
                    press(p.I) * (ug[4] - betay_avg * betaz_avg /
                                              pow2(alp_avg)); // ug[4]=uyz
    CCTK_REAL tzz = rhoenthalpyW2 * velzshift * velzshift +
                    press(p.I) * (ug[5] - betaz_avg * betaz_avg /
                                              pow2(alp_avg)); // ug[5]=uzz

    CCTK_REAL t0lowx = rhoenthalpyW2 * v_low[0] / alp_avg;
    CCTK_REAL t0lowy = rhoenthalpyW2 * v_low[1] / alp_avg;
    CCTK_REAL t0lowz = rhoenthalpyW2 * v_low[2] / alp_avg;

    /* consider magnetic field */
    {
      /* cell-centered B^i and B_i */
      const array<CCTK_REAL, 3> B_up = {Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
      const array<CCTK_REAL, 3> B_low =
          calc_vlow(B_up, gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg);
      /* b^mu */
      const CCTK_REAL bst =
          w_lorentz *
          (B_up[0] * v_low[0] + B_up[1] * v_low[1] + B_up[2] * v_low[2]) /
          alp_avg;
      const CCTK_REAL bsx =
          (B_up[0] + alp_avg * bst * w_lorentz * velxshift) / w_lorentz;
      const CCTK_REAL bsy =
          (B_up[1] + alp_avg * bst * w_lorentz * velyshift) / w_lorentz;
      const CCTK_REAL bsz =
          (B_up[2] + alp_avg * bst * w_lorentz * velzshift) / w_lorentz;
      /* b^2 */
      const CCTK_REAL bs2 = (B_up[0] * B_low[0] + B_up[1] * B_low[1] +
                             B_up[2] * B_low[2] + pow2(alp_avg * bst)) /
                            pow2(w_lorentz);

      t00 +=
          (pow2(w_lorentz / alp_avg) - 0.5 / pow2(alp_avg)) * bs2 - pow2(bst);
      t0x += (pow2(w_lorentz) * velxshift / alp_avg +
              0.5 * betax_avg / pow2(alp_avg)) *
                 bs2 -
             bst * bsx;
      t0y += (pow2(w_lorentz) * velyshift / alp_avg +
              0.5 * betay_avg / pow2(alp_avg)) *
                 bs2 -
             bst * bsy;
      t0z += (pow2(w_lorentz) * velzshift / alp_avg +
              0.5 * betaz_avg / pow2(alp_avg)) *
                 bs2 -
             bst * bsz;
      txx += (pow2(w_lorentz * velxshift) +
              0.5 * (ug[0] - pow2(betax_avg / alp_avg))) *
                 bs2 -
             pow2(bsx);
      txy += (pow2(w_lorentz) * velxshift * velyshift +
              0.5 * (ug[1] - betax_avg * betay_avg / pow2(alp_avg))) *
                 bs2 -
             bsx * bsy;
      txz += (pow2(w_lorentz) * velxshift * velzshift +
              0.5 * (ug[2] - betax_avg * betaz_avg / pow2(alp_avg))) *
                 bs2 -
             bsx * bsz;
      tyy += (pow2(w_lorentz * velyshift) +
              0.5 * (ug[3] - pow2(betay_avg / alp_avg))) *
                 bs2 -
             pow2(bsy);
      tyz += (pow2(w_lorentz) * velyshift * velzshift +
              0.5 * (ug[4] - betay_avg * betaz_avg / pow2(alp_avg))) *
                 bs2 -
             bsy * bsz;
      tzz += (pow2(w_lorentz * velzshift) +
              0.5 * (ug[5] - pow2(betaz_avg / alp_avg))) *
                 bs2 -
             pow2(bsz);

      /* beta_i */
      const array<CCTK_REAL, 3> beta_up = {betax_avg, betay_avg, betaz_avg};
      const array<CCTK_REAL, 3> beta_low = calc_vlow(
          beta_up, gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg);
      /* b_i */
      const CCTK_REAL bsx_low =
          beta_low[0] * bst + gxx_avg * bsx + gxy_avg * bsy + gxz_avg * bsz;
      const CCTK_REAL bsy_low =
          beta_low[1] * bst + gxy_avg * bsx + gyy_avg * bsy + gxz_avg * bsz;
      const CCTK_REAL bsz_low =
          beta_low[2] * bst + gxz_avg * bsx + gyz_avg * bsy + gzz_avg * bsz;
      /* T^0_i */
      t0lowx += bs2 * pow2(w_lorentz) * v_low[0] / alp_avg - bsx_low * bst;
      t0lowy += bs2 * pow2(w_lorentz) * v_low[1] / alp_avg - bsy_low * bst;
      t0lowz += bs2 * pow2(w_lorentz) * v_low[2] / alp_avg - bsz_low * bst;
    }

    /* Contract the shift with the extrinsic curvature */
    const CCTK_REAL shiftshiftk =
        betax_avg * betax_avg * kxx_avg + betay_avg * betay_avg * kyy_avg +
        betaz_avg * betaz_avg * kzz_avg +
        2.0 *
            (betax_avg * betay_avg * kxy_avg + betax_avg * betaz_avg * kxz_avg +
             betay_avg * betaz_avg * kyz_avg);
    const CCTK_REAL shiftkx =
        betax_avg * kxx_avg + betay_avg * kxy_avg + betaz_avg * kxz_avg;
    const CCTK_REAL shiftky =
        betax_avg * kxy_avg + betay_avg * kyy_avg + betaz_avg * kyz_avg;
    const CCTK_REAL shiftkz =
        betax_avg * kxz_avg + betay_avg * kyz_avg + betaz_avg * kzz_avg;

    /* Contract the matter terms with the extrinsic curvature */
    const CCTK_REAL sumTK =
        txx * kxx_avg + tyy * kyy_avg + tzz * kzz_avg +
        2.0 * (txy * kxy_avg + txz * kxz_avg + tyz * kyz_avg);

    /* Update term for tau */
    const CCTK_REAL tau_source =
        t00 * (shiftshiftk -
               (betax_avg * dx_alp + betay_avg * dy_alp + betaz_avg * dz_alp)) +
        t0x * (-dx_alp + 2.0 * shiftkx) + t0y * (-dy_alp + 2.0 * shiftky) +
        t0z * (-dz_alp + 2.0 * shiftkz) + sumTK;

    /* Contract the shift with derivatives of the metric */
    const CCTK_REAL halfshiftdgx =
        0.5 * (betax_avg * betax_avg * dx_gxx + betay_avg * betay_avg * dx_gyy +
               betaz_avg * betaz_avg * dx_gzz) +
        betax_avg * betay_avg * dx_gxy + betax_avg * betaz_avg * dx_gxz +
        betay_avg * betaz_avg * dx_gyz;

    const CCTK_REAL halfshiftdgy =
        0.5 * (betax_avg * betax_avg * dy_gxx + betay_avg * betay_avg * dy_gyy +
               betaz_avg * betaz_avg * dy_gzz) +
        betax_avg * betay_avg * dy_gxy + betax_avg * betaz_avg * dy_gxz +
        betay_avg * betaz_avg * dy_gyz;

    const CCTK_REAL halfshiftdgz =
        0.5 * (betax_avg * betax_avg * dz_gxx + betay_avg * betay_avg * dz_gyy +
               betaz_avg * betaz_avg * dz_gzz) +
        betax_avg * betay_avg * dz_gxy + betax_avg * betaz_avg * dz_gxz +
        betay_avg * betaz_avg * dz_gyz;

    /* Contract the matter with derivatives of the metric */

    const CCTK_REAL halfTdgx =
        0.5 * (txx * dx_gxx + tyy * dx_gyy + tzz * dx_gzz) + txy * dx_gxy +
        txz * dx_gxz + tyz * dx_gyz;
    const CCTK_REAL halfTdgy =
        0.5 * (txx * dy_gxx + tyy * dy_gyy + tzz * dy_gzz) + txy * dy_gxy +
        txz * dy_gxz + tyz * dy_gyz;
    const CCTK_REAL halfTdgz =
        0.5 * (txx * dz_gxx + tyy * dz_gyy + tzz * dz_gzz) + txy * dz_gxy +
        txz * dz_gxz + tyz * dz_gyz;

    const CCTK_REAL momx_source =
        t00 * (halfshiftdgx - alp_avg * dx_alp) +
        t0x * (betax_avg * dx_gxx + betay_avg * dx_gxy + betaz_avg * dx_gxz) +
        t0y * (betax_avg * dx_gxy + betay_avg * dx_gyy + betaz_avg * dx_gyz) +
        t0z * (betax_avg * dx_gxz + betay_avg * dx_gyz + betaz_avg * dx_gzz) +
        halfTdgx + t0lowx * dx_betax + t0lowy * dx_betay + t0lowz * dx_betaz;

    const CCTK_REAL momy_source =
        t00 * (halfshiftdgy - alp_avg * dy_alp) +
        t0x * (betax_avg * dy_gxx + betay_avg * dy_gxy + betaz_avg * dy_gxz) +
        t0y * (betax_avg * dy_gxy + betay_avg * dy_gyy + betaz_avg * dy_gyz) +
        t0z * (betax_avg * dy_gxz + betay_avg * dy_gyz + betaz_avg * dy_gzz) +
        halfTdgy + t0lowx * dy_betax + t0lowy * dy_betay + t0lowz * dy_betaz;

    const CCTK_REAL momz_source =
        t00 * (halfshiftdgz - alp_avg * dz_alp) +
        t0x * (betax_avg * dz_gxx + betay_avg * dz_gxy + betaz_avg * dz_gxz) +
        t0y * (betax_avg * dz_gxy + betay_avg * dz_gyy + betaz_avg * dz_gyz) +
        t0z * (betax_avg * dz_gxz + betay_avg * dz_gyz + betaz_avg * dz_gzz) +
        halfTdgz + t0lowx * dz_betax + t0lowy * dz_betay + t0lowz * dz_betaz;

    densrhs(p.I) = 0.0;
    momxrhs(p.I) = alp_avg * sqrt_detg * momx_source;
    momyrhs(p.I) = alp_avg * sqrt_detg * momy_source;
    momzrhs(p.I) = alp_avg * sqrt_detg * momz_source;
    taurhs(p.I) = alp_avg * sqrt_detg * tau_source;
  }); // end of loop over grid
}

} // namespace AsterX
