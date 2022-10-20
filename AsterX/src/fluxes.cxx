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
#include "reconstruct.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace Arith;

// Calculate the fluxes in direction `dir`. This function is more
// complex because it has to handle any direction, but as reward,
// there is only one function, not three.
template <int dir> void CalcFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};

  // Get the grid function pointers for fluxes in direction `dir`
  const array<GF3D2<CCTK_REAL>, dim> fluxdenss = {fxdens, fydens, fzdens};
  const array<GF3D2<CCTK_REAL>, dim> fluxmomxs = {fxmomx, fymomx, fzmomx};
  const array<GF3D2<CCTK_REAL>, dim> fluxmomys = {fxmomy, fymomy, fzmomy};
  const array<GF3D2<CCTK_REAL>, dim> fluxmomzs = {fxmomz, fymomz, fzmomz};
  const array<GF3D2<CCTK_REAL>, dim> fluxtaus = {fxtau, fytau, fztau};
  const array<GF3D2<CCTK_REAL>, dim> fluxBxs = {fxBx, fyBx, fzBx};
  const array<GF3D2<CCTK_REAL>, dim> fluxBys = {fxBy, fyBy, fzBy};
  const array<GF3D2<CCTK_REAL>, dim> fluxBzs = {fxBz, fyBz, fzBz};

  // fdens^i = rho vel^i
  // fmom^i_j = mom_j vel^i + delta^i_j press
  // ftau^i = (tau + press) vel^i

  reconstruction_t reconstruction;
  if (CCTK_EQUALS(reconstruction_method, "Godunov"))
    reconstruction = reconstruction_t::Godunov;
  else if (CCTK_EQUALS(reconstruction_method, "minmod"))
    reconstruction = reconstruction_t::minmod;
  else if (CCTK_EQUALS(reconstruction_method, "monocentral"))
    reconstruction = reconstruction_t::monocentral;
  else if (CCTK_EQUALS(reconstruction_method, "ppm"))
    reconstruction = reconstruction_t::ppm;
  else
    CCTK_ERROR("Unknown value for parameter \"reconstruction_method\"");

  switch (reconstruction) {
  case reconstruction_t::Godunov:
    assert(cctk_nghostzones[dir] >= 1);
    break;
  case reconstruction_t::minmod:
    assert(cctk_nghostzones[dir] >= 2);
    break;
  case reconstruction_t::monocentral:
    assert(cctk_nghostzones[dir] >= 2);
    break;
  case reconstruction_t::ppm:
    assert(cctk_nghostzones[dir] >= 3);
    break;
  }

  const auto reconstruct_pt =
      [=] CCTK_DEVICE(const GF3D2<const CCTK_REAL> &var, const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            return reconstruct(var, p, reconstruction, dir);
          };

  const auto eigenvalues =
      [=] CCTK_DEVICE(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                      array<CCTK_REAL, 2> vel, array<CCTK_REAL, 2> rho,
                      array<CCTK_REAL, 2> cs2, array<CCTK_REAL, 2> w_lor,
                      array<CCTK_REAL, 2> h,
                      array<CCTK_REAL, 2> bsq) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // computing characteristics for the minus side
        // See Eq. (28) of Giacomazzo & Rezzolla (2007) with b^i=0
        array<CCTK_REAL, 3> a_m = {
            (bsq[0] + cs2[0] * h[0] * rho[0]) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2[0]) * h[0] * rho[0] *
                    pow2(beta_avg - alp_avg * vel[0]) * pow2(w_lor[0]),

            2 * beta_avg * (bsq[0] + cs2[0] * h[0] * rho[0]) -
                2 * (-1 + cs2[0]) * h[0] * rho[0] *
                    (beta_avg - alp_avg * vel[0]) * pow2(w_lor[0]),

            bsq[0] + h[0] * rho[0] *
                         (cs2[0] + pow2(w_lor[0]) - cs2[0] * pow2(w_lor[0]))};

        CCTK_REAL det_m = pow2(a_m[1]) - 4.0 * a_m[2] * a_m[0];
        if (det_m < 0.0)
          det_m = 0.0;

        array<CCTK_REAL, 4> lambda_m = {
            ((-a_m[1] + sqrt(det_m)) / (2.0 * a_m[2])) / alp_avg,
            ((-a_m[1] + sqrt(det_m)) / (2.0 * a_m[2])) / alp_avg,
            ((-a_m[1] - sqrt(det_m)) / (2.0 * a_m[2])) / alp_avg,
            ((-a_m[1] - sqrt(det_m)) / (2.0 * a_m[2])) / alp_avg};

        // computing characteristics for the plus side

        array<CCTK_REAL, 3> a_p = {
            (bsq[1] + cs2[1] * h[1] * rho[1]) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2[1]) * h[1] * rho[1] *
                    pow2(beta_avg - alp_avg * vel[1]) * pow2(w_lor[1]),

            2 * beta_avg * (bsq[1] + cs2[1] * h[1] * rho[1]) -
                2 * (-1 + cs2[1]) * h[1] * rho[1] *
                    (beta_avg - alp_avg * vel[1]) * pow2(w_lor[1]),

            bsq[1] + h[1] * rho[1] *
                         (cs2[1] + pow2(w_lor[1]) - cs2[1] * pow2(w_lor[1]))};

        CCTK_REAL det_p = pow2(a_p[1]) - 4.0 * a_p[2] * a_p[0];
        if (det_p < 0.0)
          det_p = 0.0;

        array<CCTK_REAL, 4> lambda_p = {
            ((-a_p[1] + sqrt(det_p)) / (2.0 * a_p[2])) / alp_avg,
            ((-a_p[1] + sqrt(det_p)) / (2.0 * a_p[2])) / alp_avg,
            ((-a_p[1] - sqrt(det_p)) / (2.0 * a_p[2])) / alp_avg,
            ((-a_p[1] - sqrt(det_p)) / (2.0 * a_p[2])) / alp_avg};

        // 2D array containing characteristics for left (minus) and right (plus)
        // sides
        array<array<CCTK_REAL, 4>, 2> lambda = {lambda_m, lambda_p};
        return lambda;
      };

  const auto calcflux =
      [=] CCTK_DEVICE(array<array<CCTK_REAL, 4>, 2> lam,
                      array<CCTK_REAL, 2> var, array<CCTK_REAL, 2> flux)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
            const CCTK_REAL charmax =
                max({0.0, fabs(lam[0][0]), fabs(lam[0][1]), fabs(lam[0][2]),
                     fabs(lam[0][3]), fabs(lam[1][0]), fabs(lam[1][1]),
                     fabs(lam[1][2]), fabs(lam[1][3])});

            CCTK_REAL llf =
                0.5 * ((flux[0] + flux[1]) - charmax * (var[1] - var[0]));
            // return dA * llf;
            return llf;
          };

  grid.loop_int_device<
      face_centred[0], face_centred[1],
      face_centred
          [2]>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    // Reconstruct primitives from the cells on left (indice 0) and right
    // (indice 1) side of this face rc = reconstructed variables or computed
    // from reconstructed variables
    const array<CCTK_REAL, 2> rho_rc = reconstruct_pt(rho, p);
    const array<CCTK_REAL, 2> velx_rc = reconstruct_pt(velx, p);
    const array<CCTK_REAL, 2> vely_rc = reconstruct_pt(vely, p);
    const array<CCTK_REAL, 2> velz_rc = reconstruct_pt(velz, p);
    const array<CCTK_REAL, 2> eps_rc = reconstruct_pt(eps, p);
    const array<CCTK_REAL, 2> Bx_rc = reconstruct_pt(Bvecx, p);
    const array<CCTK_REAL, 2> By_rc = reconstruct_pt(Bvecy, p);
    const array<CCTK_REAL, 2> Bz_rc = reconstruct_pt(Bvecz, p);

    const array<array<CCTK_REAL, 2>, 3> vels_rc = {velx_rc, vely_rc, velz_rc};
    const array<CCTK_REAL, 2> vel_rc = vels_rc[dir];

    const array<array<CCTK_REAL, 2>, 3> Bs_rc = {Bx_rc, By_rc, Bz_rc};
    const array<CCTK_REAL, 2> B_rc = Bs_rc[dir];

    // TODO: to reconstruct w_lorentz*vel or 4-velocity u_i

    // Computing metric components
    CCTK_REAL alp_avg = calc_avg_v2f(alp, p, dir);
    CCTK_REAL betax_avg = calc_avg_v2f(betax, p, dir);
    CCTK_REAL betay_avg = calc_avg_v2f(betay, p, dir);
    CCTK_REAL betaz_avg = calc_avg_v2f(betaz, p, dir);
    CCTK_REAL gxx_avg = calc_avg_v2f(gxx, p, dir);
    CCTK_REAL gxy_avg = calc_avg_v2f(gxy, p, dir);
    CCTK_REAL gxz_avg = calc_avg_v2f(gxz, p, dir);
    CCTK_REAL gyy_avg = calc_avg_v2f(gyy, p, dir);
    CCTK_REAL gyz_avg = calc_avg_v2f(gyz, p, dir);
    CCTK_REAL gzz_avg = calc_avg_v2f(gzz, p, dir);

    const array<CCTK_REAL, 3> betas_avg = {betax_avg, betay_avg, betaz_avg};
    const CCTK_REAL beta_avg = betas_avg[dir];

    // TODO: Compute pressure based on user-specified EOS.
    // Currently, computing press for classical ideal gas from reconstructed
    // vars

    const array<CCTK_REAL, 2> press_rc = {eps_rc[0] * rho_rc[0] * (gamma - 1),
                                          eps_rc[1] * rho_rc[1] * (gamma - 1)};

    // Determinant of spatial metric
    const smat<CCTK_REAL, 3> g_avg{gxx_avg, gxy_avg, gxz_avg,
                                           gyy_avg, gyz_avg, gzz_avg};
    const CCTK_REAL detg_avg = calc_det(g_avg);
    const CCTK_REAL sqrt_detg = sqrt(detg_avg);

    // Upper metric
    const smat<CCTK_REAL, 3> ug_avg = calc_inv(g_avg, detg_avg);

    // Variable for either uxx, uyy or uzz depending on the direction
    const CCTK_REAL u_avg = ug_avg(dir, dir);

    // v_j
    const array<CCTK_REAL, 2> vlowx_rc = {
        gxx_avg * velx_rc[0] + gxy_avg * vely_rc[0] + gxz_avg * velz_rc[0],
        gxx_avg * velx_rc[1] + gxy_avg * vely_rc[1] + gxz_avg * velz_rc[1]};

    const array<CCTK_REAL, 2> vlowy_rc = {
        gxy_avg * velx_rc[0] + gyy_avg * vely_rc[0] + gyz_avg * velz_rc[0],
        gxy_avg * velx_rc[1] + gyy_avg * vely_rc[1] + gyz_avg * velz_rc[1]};

    const array<CCTK_REAL, 2> vlowz_rc = {
        gxz_avg * velx_rc[0] + gyz_avg * vely_rc[0] + gzz_avg * velz_rc[0],
        gxz_avg * velx_rc[1] + gyz_avg * vely_rc[1] + gzz_avg * velz_rc[1]};

    // Computing the contravariant coordinate velocity
    // vtilde^i = alpha*v^i - beta^i using the reconstructed variables
    const array<CCTK_REAL, 2> vtildex_rc = {alp_avg * velx_rc[0] - betax_avg,
                                            alp_avg * velx_rc[1] - betax_avg};

    const array<CCTK_REAL, 2> vtildey_rc = {alp_avg * vely_rc[0] - betay_avg,
                                            alp_avg * vely_rc[1] - betay_avg};

    const array<CCTK_REAL, 2> vtildez_rc = {alp_avg * velz_rc[0] - betaz_avg,
                                            alp_avg * velz_rc[1] - betaz_avg};

    const array<array<CCTK_REAL, 2>, 3> vtildes_rc = {vtildex_rc, vtildey_rc,
                                                      vtildez_rc};

    const array<CCTK_REAL, 2> vtilde_rc = vtildes_rc[dir];

    // FIXME: consider using the function from utils.hxx
    // Computing w_lorentz using reconstructed variables
    const array<CCTK_REAL, 2> w_lorentz_rc = {
        1.0 / sqrt(1 - (vlowx_rc[0] * velx_rc[0] + vlowy_rc[0] * vely_rc[0] +
                        vlowz_rc[0] * velz_rc[0])),
        1.0 / sqrt(1 - (vlowx_rc[1] * velx_rc[1] + vlowy_rc[1] * vely_rc[1] +
                        vlowz_rc[1] * velz_rc[1]))};

    // Computing cs2 for ideal gas EOS using reconstructed variables
    const array<CCTK_REAL, 2> cs2_rc = {
        (gamma - 1.0) * eps_rc[0] / (eps_rc[0] + 1.0 / gamma),
        (gamma - 1.0) * eps_rc[1] / (eps_rc[1] + 1.0 / gamma)};

    // Computing enthalpy h for ideal gas EOS using reconstructed variables
    const array<CCTK_REAL, 2> h_rc = {1.0 + eps_rc[0] + press_rc[0] / rho_rc[0],
                                      1.0 + eps_rc[1] +
                                          press_rc[1] / rho_rc[1]};

    // Computing the covariant magnetic field measured by the Eulerian observer
    // using the reconstructed variables
    const array<CCTK_REAL, 2> Blowx_rc = {
        gxx_avg * Bx_rc[0] + gxy_avg * By_rc[0] + gxz_avg * Bz_rc[0],
        gxx_avg * Bx_rc[1] + gxy_avg * By_rc[1] + gxz_avg * Bz_rc[1]};

    const array<CCTK_REAL, 2> Blowy_rc = {
        gxy_avg * Bx_rc[0] + gyy_avg * By_rc[0] + gyz_avg * Bz_rc[0],
        gxy_avg * Bx_rc[1] + gyy_avg * By_rc[1] + gyz_avg * Bz_rc[1]};

    const array<CCTK_REAL, 2> Blowz_rc = {
        gxz_avg * Bx_rc[0] + gyz_avg * By_rc[0] + gzz_avg * Bz_rc[0],
        gxz_avg * Bx_rc[1] + gyz_avg * By_rc[1] + gzz_avg * Bz_rc[1]};

    const array<CCTK_REAL, 2> B2_rc = {
        Bx_rc[0] * Blowx_rc[0] + By_rc[0] * Blowy_rc[0] +
            Bz_rc[0] * Blowz_rc[0],
        Bx_rc[1] * Blowx_rc[1] + By_rc[1] * Blowy_rc[1] +
            Bz_rc[1] * Blowz_rc[1]};

    // Computing the magnetic field measured by the observer comoving with the
    // fluid using the reconstructed variables
    const array<CCTK_REAL, 2> alpha_b0_rc = {
        w_lorentz_rc[0] * (Bx_rc[0] * vlowx_rc[0] + By_rc[0] * vlowy_rc[0] +
                           Bz_rc[0] * vlowz_rc[0]),
        w_lorentz_rc[1] * (Bx_rc[1] * vlowx_rc[1] + By_rc[1] * vlowy_rc[1] +
                           Bz_rc[1] * vlowz_rc[1])};

    const array<CCTK_REAL, 2> blowx_rc = {
        Blowx_rc[0] / w_lorentz_rc[0] + alpha_b0_rc[0] * vlowx_rc[0],
        Blowx_rc[1] / w_lorentz_rc[1] + alpha_b0_rc[1] * vlowx_rc[1]};

    const array<CCTK_REAL, 2> blowy_rc = {
        Blowy_rc[0] / w_lorentz_rc[0] + alpha_b0_rc[0] * vlowy_rc[0],
        Blowy_rc[1] / w_lorentz_rc[1] + alpha_b0_rc[1] * vlowy_rc[1]};

    const array<CCTK_REAL, 2> blowz_rc = {
        Blowz_rc[0] / w_lorentz_rc[0] + alpha_b0_rc[0] * vlowz_rc[0],
        Blowz_rc[1] / w_lorentz_rc[1] + alpha_b0_rc[1] * vlowz_rc[1]};

    const array<CCTK_REAL, 2> bsq_rc = {
        (B2_rc[0] + pow2(alpha_b0_rc[0])) / pow2(w_lorentz_rc[0]),
        (B2_rc[1] + pow2(alpha_b0_rc[1])) / pow2(w_lorentz_rc[1])};

    // Auxiliary variables to compute the conservative variables and their
    // fluxes
    const array<CCTK_REAL, 2> sqrt_detg_press_plus_pmag_rc = {
        // sqrt(g)*(p_gas + p_mag)
        sqrt_detg * (press_rc[0] + 0.5 * bsq_rc[0]),
        sqrt_detg * (press_rc[1] + 0.5 * bsq_rc[1])};

    const array<CCTK_REAL, 2> alp_sqrt_detg_press_plus_pmag_rc = {
        // alp*sqrt(g)*(p_gas + p_mag)
        alp_avg * sqrt_detg_press_plus_pmag_rc[0],
        alp_avg * sqrt_detg_press_plus_pmag_rc[1]};

    const array<CCTK_REAL, 2> alp_sqrt_detg_B_over_w_lorentz_rc = {
        // alp*sqrt(g)*B^i / W
        alp_avg * sqrt_detg * B_rc[0] / w_lorentz_rc[0],
        alp_avg * sqrt_detg * B_rc[1] / w_lorentz_rc[1]};

    // Computing conservatives from primitives
    const array<CCTK_REAL, 2> dens_rc = {
        // sqrt(g)*rho*W
        sqrt_detg * rho_rc[0] * w_lorentz_rc[0],
        sqrt_detg * rho_rc[1] * w_lorentz_rc[1]};

    const array<CCTK_REAL, 2> dens_h_W_rc = {
        // D*h*W = sqrt(g)*rho*h*W^2
        dens_rc[0] * h_rc[0] * w_lorentz_rc[0],
        dens_rc[1] * h_rc[1] * w_lorentz_rc[1]};

    const array<CCTK_REAL, 2> dens_h_W_plus_sqrt_detg_W2b2_rc = {
        // sqrt(g)*( rho*h*W^2 + b^2*W^2 ) = sqrt(g)( rho*h*W^2 + (alp^2*b0^2) +
        // B^2 )
        dens_h_W_rc[0] + sqrt_detg * (pow2(alpha_b0_rc[0]) + B2_rc[0]),
        dens_h_W_rc[1] + sqrt_detg * (pow2(alpha_b0_rc[1]) + B2_rc[1])};

    const array<CCTK_REAL, 2> momx_rc = {
        dens_h_W_plus_sqrt_detg_W2b2_rc[0] * vlowx_rc[0] -
            sqrt_detg * alpha_b0_rc[0] * blowx_rc[0],
        dens_h_W_plus_sqrt_detg_W2b2_rc[1] * vlowx_rc[1] -
            sqrt_detg * alpha_b0_rc[1] * blowx_rc[1]};

    const array<CCTK_REAL, 2> momy_rc = {
        dens_h_W_plus_sqrt_detg_W2b2_rc[0] * vlowy_rc[0] -
            sqrt_detg * alpha_b0_rc[0] * blowy_rc[0],
        dens_h_W_plus_sqrt_detg_W2b2_rc[1] * vlowy_rc[1] -
            sqrt_detg * alpha_b0_rc[1] * blowy_rc[1]};

    const array<CCTK_REAL, 2> momz_rc = {
        dens_h_W_plus_sqrt_detg_W2b2_rc[0] * vlowz_rc[0] -
            sqrt_detg * alpha_b0_rc[0] * blowz_rc[0],
        dens_h_W_plus_sqrt_detg_W2b2_rc[1] * vlowz_rc[1] -
            sqrt_detg * alpha_b0_rc[1] * blowz_rc[1]};

    // FIXME: B^2 = W^2·b^2 - (alpha·b^0)^2, is that true?
    const array<CCTK_REAL, 2> tau_rc = {
        dens_h_W_rc[0] - dens_rc[0] - sqrt_detg_press_plus_pmag_rc[0] +
            sqrt_detg * B2_rc[0],
        dens_h_W_rc[1] - dens_rc[1] - sqrt_detg_press_plus_pmag_rc[1] +
            sqrt_detg * B2_rc[1]};

    const array<CCTK_REAL, 2> Btildex_rc = {sqrt_detg * Bx_rc[0],
                                            sqrt_detg * Bx_rc[1]};

    const array<CCTK_REAL, 2> Btildey_rc = {sqrt_detg * By_rc[0],
                                            sqrt_detg * By_rc[1]};

    const array<CCTK_REAL, 2> Btildez_rc = {sqrt_detg * Bz_rc[0],
                                            sqrt_detg * Bz_rc[1]};

    // Computing fluxes of conserved variables
    const array<CCTK_REAL, 2> flux_dens = {dens_rc[0] * vtilde_rc[0],
                                           dens_rc[1] * vtilde_rc[1]};

    const array<CCTK_REAL, 2> flux_momx = {
        momx_rc[0] * vtilde_rc[0] +
            (dir == 0) * alp_sqrt_detg_press_plus_pmag_rc[0] -
            alp_sqrt_detg_B_over_w_lorentz_rc[0] * blowx_rc[0],
        momx_rc[1] * vtilde_rc[1] +
            (dir == 0) * alp_sqrt_detg_press_plus_pmag_rc[1] -
            alp_sqrt_detg_B_over_w_lorentz_rc[1] * blowx_rc[1]};

    const array<CCTK_REAL, 2> flux_momy = {
        momy_rc[0] * vtilde_rc[0] +
            (dir == 1) * alp_sqrt_detg_press_plus_pmag_rc[0] -
            alp_sqrt_detg_B_over_w_lorentz_rc[0] * blowy_rc[0],
        momy_rc[1] * vtilde_rc[1] +
            (dir == 1) * alp_sqrt_detg_press_plus_pmag_rc[1] -
            alp_sqrt_detg_B_over_w_lorentz_rc[1] * blowy_rc[1]};

    const array<CCTK_REAL, 2> flux_momz = {
        momz_rc[0] * vtilde_rc[0] +
            (dir == 2) * alp_sqrt_detg_press_plus_pmag_rc[0] -
            alp_sqrt_detg_B_over_w_lorentz_rc[0] * blowz_rc[0],
        momz_rc[1] * vtilde_rc[1] +
            (dir == 2) * alp_sqrt_detg_press_plus_pmag_rc[1] -
            alp_sqrt_detg_B_over_w_lorentz_rc[1] * blowz_rc[1]};

    const array<CCTK_REAL, 2> flux_tau = {
        tau_rc[0] * vtilde_rc[0] +
            alp_sqrt_detg_press_plus_pmag_rc[0] * vel_rc[0] -
            alpha_b0_rc[0] * alp_sqrt_detg_B_over_w_lorentz_rc[0],
        tau_rc[1] * vtilde_rc[1] +
            alp_sqrt_detg_press_plus_pmag_rc[1] * vel_rc[1] -
            alpha_b0_rc[1] * alp_sqrt_detg_B_over_w_lorentz_rc[1]};

    const array<CCTK_REAL, 2> flux_Btildex = {
        // (0, Ez, -Ey)
        (dir ==
         1) * (Btildex_rc[0] * vtildey_rc[0] - Btildey_rc[0] * vtildex_rc[0]) -
            (dir == 2) *
                (Btildez_rc[0] * vtildex_rc[0] - Btildex_rc[0] * vtildez_rc[0]),
        (dir ==
         1) * (Btildex_rc[1] * vtildey_rc[1] - Btildey_rc[1] * vtildex_rc[1]) -
            (dir == 2) * (Btildez_rc[1] * vtildex_rc[1] -
                          Btildex_rc[1] * vtildez_rc[1])};

    const array<CCTK_REAL, 2> flux_Btildey = {
        // (-Ez, 0, Ex)
        -(dir ==
          0) * (Btildex_rc[0] * vtildey_rc[0] - Btildey_rc[0] * vtildex_rc[0]) +
            (dir == 2) *
                (Btildey_rc[0] * vtildez_rc[0] - Btildez_rc[0] * vtildey_rc[0]),
        -(dir ==
          0) * (Btildex_rc[1] * vtildey_rc[1] - Btildey_rc[1] * vtildex_rc[1]) +
            (dir == 2) * (Btildey_rc[1] * vtildez_rc[1] -
                          Btildez_rc[1] * vtildey_rc[1])};

    const array<CCTK_REAL, 2> flux_Btildez = {
        // (Ey, -Ex, 0)
        (dir ==
         0) * (Btildez_rc[0] * vtildex_rc[0] - Btildex_rc[0] * vtildez_rc[0]) -
            (dir == 1) *
                (Btildey_rc[0] * vtildez_rc[0] - Btildez_rc[0] * vtildey_rc[0]),
        (dir ==
         0) * (Btildez_rc[1] * vtildex_rc[1] - Btildex_rc[1] * vtildez_rc[1]) -
            (dir == 1) * (Btildey_rc[1] * vtildez_rc[1] -
                          Btildez_rc[1] * vtildey_rc[1])};

    array<array<CCTK_REAL, 4>, 2> lambda =
        eigenvalues(alp_avg, beta_avg, u_avg, vel_rc, rho_rc, cs2_rc,
                    w_lorentz_rc, h_rc, bsq_rc);

    fluxdenss[dir](p.I) = calcflux(lambda, dens_rc, flux_dens);
    fluxmomxs[dir](p.I) = calcflux(lambda, momx_rc, flux_momx);
    fluxmomys[dir](p.I) = calcflux(lambda, momy_rc, flux_momy);
    fluxmomzs[dir](p.I) = calcflux(lambda, momz_rc, flux_momz);
    fluxtaus[dir](p.I) = calcflux(lambda, tau_rc, flux_tau);
    fluxBxs[dir](p.I) = (dir != 0) * calcflux(lambda, Btildex_rc, flux_Btildex);
    fluxBys[dir](p.I) = (dir != 1) * calcflux(lambda, Btildey_rc, flux_Btildey);
    fluxBzs[dir](p.I) = (dir != 2) * calcflux(lambda, Btildez_rc, flux_Btildez);
  });
}

void CalcAuxForAvecPsi(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* interpolate A to vertices */
        const CCTK_REAL Ax_vert = calc_avg_e2v(Avec_x, p, 0);
        const CCTK_REAL Ay_vert = calc_avg_e2v(Avec_y, p, 1);
        const CCTK_REAL Az_vert = calc_avg_e2v(Avec_z, p, 2);

        const smat<CCTK_REAL, 3> g{gxx(p.I), gxy(p.I), gxz(p.I),
                                           gyy(p.I), gyz(p.I), gzz(p.I)};
        const CCTK_REAL detg = calc_det(g);
        const smat<CCTK_REAL, 3> ug = calc_inv(g, detg);

        const CCTK_REAL Axup =
            ug(0,0) * Ax_vert + ug(0,1) * Ay_vert + ug(0,2) * Az_vert;
        const CCTK_REAL Ayup =
            ug(0,1) * Ax_vert + ug(1,1) * Ay_vert + ug(1,2) * Az_vert;
        const CCTK_REAL Azup =
            ug(0,2) * Ax_vert + ug(1,2) * Ay_vert + ug(2,2) * Az_vert;

        const CCTK_REAL beta_Avec =
            betax(p.I) * Ax_vert + betay(p.I) * Ay_vert + betaz(p.I) * Az_vert;

        Fx(p.I) = alp(p.I) * sqrt(detg) * Axup;
        Fy(p.I) = alp(p.I) * sqrt(detg) * Ayup;
        Fz(p.I) = alp(p.I) * sqrt(detg) * Azup;
        Fbetax(p.I) = betax(p.I) * Psi(p.I);
        Fbetay(p.I) = betay(p.I) * Psi(p.I);
        Fbetaz(p.I) = betaz(p.I) * Psi(p.I);
        G(p.I) = alp(p.I) * Psi(p.I) / sqrt(detg) - beta_Avec;
      });
}

extern "C" void AsterX_Fluxes(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  CalcFlux<0>(cctkGH);
  CalcFlux<1>(cctkGH);
  CalcFlux<2>(cctkGH);

  /* Set auxiliary variables for the rhs of A and Psi  */
  CalcAuxForAvecPsi(cctkGH);
}

} // namespace AsterX
