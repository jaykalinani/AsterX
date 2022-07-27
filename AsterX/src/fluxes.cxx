#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "utils.hxx"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace AsterX {
using namespace std;
using namespace Loop;

enum class reconstruction_t { Godunov, minmod, monocentral, ppm };

namespace {
template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T minmod(const T &x,
                                                                   const T &y) {
  if (signbit(x) != signbit(y))
    return T(0);
  if (fabs(x) < fabs(y))
    return x;
  else
    return y;
}

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
monocentral(const T &x, const T &y) {
  if (sgn(x) != sgn(y))
    return 0;
  else
    return sgn(x) * min(2.0 * fabs(x), min(2.0 * fabs(y), fabs(x + y) / 2));
}
} // namespace

// Calculate the fluxes in direction `dir`. This function is more
// complex because it has to handle any direction, but as reward,
// there is only one function, not three.
template <int dir> void CalcFlux(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Fluxes;
  DECLARE_CCTK_PARAMETERS;

  static_assert(dir >= 0 && dir < 3, "");

  // const array<CCTK_REAL, dim> dx = {CCTK_DELTA_SPACE(0), CCTK_DELTA_SPACE(1),
  //                                   CCTK_DELTA_SPACE(2)};
  // const CCTK_REAL dV = dx[0] * dx[1] * dx[2]; // cell volume
  // const CCTK_REAL dA = dV / dx[dir];          // face area

  // Cell centred grid functions
  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  constexpr array<int, dim> vertex_centred = {0, 0, 0};
  const GF3D2layout gf_layout_cell(cctkGH, cell_centred);
  const GF3D2layout gf_layout_vertex(cctkGH, vertex_centred);

  const GF3D2<const CCTK_REAL> gf_gxx(gf_layout_vertex, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(gf_layout_vertex, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(gf_layout_vertex, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(gf_layout_vertex, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(gf_layout_vertex, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(gf_layout_vertex, gzz);

  const GF3D2<const CCTK_REAL> gf_alp(gf_layout_vertex, alp);
  const GF3D2<const CCTK_REAL> gf_betax(gf_layout_vertex, betax);
  const GF3D2<const CCTK_REAL> gf_betay(gf_layout_vertex, betay);
  const GF3D2<const CCTK_REAL> gf_betaz(gf_layout_vertex, betaz);

  const GF3D2<const CCTK_REAL> gf_dens(gf_layout_cell, dens);
  const GF3D2<const CCTK_REAL> gf_momx(gf_layout_cell, momx);
  const GF3D2<const CCTK_REAL> gf_momy(gf_layout_cell, momy);
  const GF3D2<const CCTK_REAL> gf_momz(gf_layout_cell, momz);
  const GF3D2<const CCTK_REAL> gf_tau(gf_layout_cell, tau);

  const GF3D2<const CCTK_REAL> gf_rho(gf_layout_cell, rho);
  const GF3D2<const CCTK_REAL> gf_velx(gf_layout_cell, velx);
  const GF3D2<const CCTK_REAL> gf_vely(gf_layout_cell, vely);
  const GF3D2<const CCTK_REAL> gf_velz(gf_layout_cell, velz);
  const GF3D2<const CCTK_REAL> gf_press(gf_layout_cell, press);
  const GF3D2<const CCTK_REAL> gf_eps(gf_layout_cell, eps);

  // FIXME: are Bvecx, Bvecy, Bvecz really GFs? I can't find where they are defined
  const GF3D2<const CCTK_REAL> gf_Bx(gf_layout_cell, Bvecx);
  const GF3D2<const CCTK_REAL> gf_By(gf_layout_cell, Bvecy);
  const GF3D2<const CCTK_REAL> gf_Bz(gf_layout_cell, Bvecz);

  // Face-centred grid functions (in direction `dir`)
  constexpr array<int, dim> face_centred = {!(dir == 0), !(dir == 1),
                                            !(dir == 2)};
  const GF3D2layout gf_fluxlayout(cctkGH, face_centred);

  // Get the grid function pointers for fluxes in direction `dir`
  const array<CCTK_REAL *, dim> fluxdenss = {fxdens, fydens, fzdens};
  const array<CCTK_REAL *, dim> fluxmomxs = {fxmomx, fymomx, fzmomx};
  const array<CCTK_REAL *, dim> fluxmomys = {fxmomy, fymomy, fzmomy};
  const array<CCTK_REAL *, dim> fluxmomzs = {fxmomz, fymomz, fzmomz};
  const array<CCTK_REAL *, dim> fluxtaus = {fxtau, fytau, fztau};
  const GF3D2<CCTK_REAL> gf_fluxdens(gf_fluxlayout, fluxdenss[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomx(gf_fluxlayout, fluxmomxs[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomy(gf_fluxlayout, fluxmomys[dir]);
  const GF3D2<CCTK_REAL> gf_fluxmomz(gf_fluxlayout, fluxmomzs[dir]);
  const GF3D2<CCTK_REAL> gf_fluxtau(gf_fluxlayout, fluxtaus[dir]);

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

  constexpr auto DI = PointDesc::DI;

  const auto reconstruct =
      [=] CCTK_DEVICE(const GF3D2<const CCTK_REAL> &gf_var,
                      const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Neighbouring "plus" and "minus" cell indices
        const auto Immm = p.I - 3 * DI[dir];
        const auto Imm = p.I - 2 * DI[dir];
        const auto Im = p.I - DI[dir];
        const auto Ip = p.I;
        const auto Ipp = p.I + DI[dir];
        const auto Ippp = p.I + 2 * DI[dir];

        switch (reconstruction) {
        case reconstruction_t::Godunov: {
          CCTK_REAL var_m = gf_var(Im);
          CCTK_REAL var_p = gf_var(Ip);
          return array<CCTK_REAL, 2>{var_m, var_p};
        }

        case reconstruction_t::minmod: {
          // reconstructs values of Im and Ip at the common face between these
          // two cells
          CCTK_REAL var_slope_p = gf_var(Ipp) - gf_var(Ip);
          CCTK_REAL var_slope_c = gf_var(Ip) - gf_var(Im);
          CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
          // reconstructed Im on its "plus/right" side
          CCTK_REAL var_m = gf_var(Im) + minmod(var_slope_c, var_slope_m) / 2;
          // reconstructed Ip on its "minus/left" side
          CCTK_REAL var_p = gf_var(Ip) - minmod(var_slope_p, var_slope_c) / 2;
          return array<CCTK_REAL, 2>{var_m, var_p};
        }

        case reconstruction_t::monocentral: {
          // reconstructs values of Im and Ip at the common face between these
          // two cells
          // reconstructed Im on its "plus/right" side
          CCTK_REAL var_slope_p = gf_var(Ip) - gf_var(Im);
          CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
          CCTK_REAL var_m =
              gf_var(Im) + monocentral(var_slope_p, var_slope_m) / 2;
          // reconstructed Ip on its "minus/left" side
          var_slope_p = gf_var(Ipp) - gf_var(Ip);
          var_slope_m = gf_var(Ip) - gf_var(Im);
          CCTK_REAL var_p =
              gf_var(Ip) - monocentral(var_slope_p, var_slope_m) / 2;
          return array<CCTK_REAL, 2>{var_m, var_p};
        }

        case reconstruction_t::ppm: {
          // Usually recon methods return the left and right states of a given
          // cell face, ie (i-1/2-eps) and (i-1/2+eps).
          // But ppm return the states at the left and right
          // faces of a given cell, ie (i-1/2) and (i+1/2).
          // So, to get left and right states from (i-1/2),
          // we first apply ppm to cell i and keep left face (i-1/2+eps)
          // and then apply ppm to cell i-1 and keep right face (i-1/2-eps)

          // Start calculating left (i-1/2) and right (i+1/2) faces unique
          // values: Eq. (A1) in https://arxiv.org/pdf/astro-ph/0503420.pdf with
          // 1/8 --> 1/6 Equiv. to Eq. (1.6) in C&W (1984)
          CCTK_REAL var_slope_p = gf_var(Ip) - gf_var(Im);
          CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
          CCTK_REAL grad_m = monocentral(var_slope_p, var_slope_m);
          var_slope_p = gf_var(Ipp) - gf_var(Ip);
          var_slope_m = gf_var(Ip) - gf_var(Im);
          CCTK_REAL grad_p = monocentral(var_slope_p, var_slope_m);
          CCTK_REAL left_face =
              (gf_var(Im) + gf_var(Ip)) / 2.0 + (grad_m - grad_p) / 6.0;

          var_slope_p = gf_var(Ipp) - gf_var(Ip);
          var_slope_m = gf_var(Ip) - gf_var(Im);
          grad_m = monocentral(var_slope_p, var_slope_m);
          var_slope_p = gf_var(Ippp) - gf_var(Ipp);
          var_slope_m = gf_var(Ipp) - gf_var(Ip);
          grad_p = monocentral(var_slope_p, var_slope_m);
          CCTK_REAL right_face =
              (gf_var(Ip) + gf_var(Ipp)) / 2.0 + (grad_m - grad_p) / 6.0;

          // Now apply conditions in Eq. (1.11) of C&W (1984)
          CCTK_REAL qa = (right_face - gf_var(Ip)) * (gf_var(Ip) - left_face);
          CCTK_REAL qd = (right_face - left_face);
          CCTK_REAL qe = 6.0 * (gf_var(Ip) - (left_face + right_face) / 2.0);
          if (qa <= 0.) {
            left_face = gf_var(Ip);
            right_face = gf_var(Ip);
          }
          if (qd * (qd - qe) < 0.0)
            left_face = 3.0 * gf_var(Ip) - 2.0 * right_face;
          if (qd * (qd + qe) < 0.0)
            right_face = 3.0 * gf_var(Ip) - 2.0 * left_face;

          // Keep left value of cell i as i-1/2+eps
          CCTK_REAL var_p = left_face;

          // Start calculating left (i-1-1/2) and right (i-1+1/2) faces unique
          // values: Eq. (A1) in https://arxiv.org/pdf/astro-ph/0503420.pdf with
          // 1/8 --> 1/6 Equiv. to Eq. (1.6) in C&W (1984)
          var_slope_p = gf_var(Im) - gf_var(Imm);
          var_slope_m = gf_var(Imm) - gf_var(Immm);
          grad_m = monocentral(var_slope_p, var_slope_m);
          var_slope_p = gf_var(Ip) - gf_var(Im);
          var_slope_m = gf_var(Im) - gf_var(Imm);
          grad_p = monocentral(var_slope_p, var_slope_m);
          left_face =
              (gf_var(Imm) + gf_var(Im)) / 2.0 + (grad_m - grad_p) / 6.0;

          var_slope_p = gf_var(Ip) - gf_var(Im);
          var_slope_m = gf_var(Im) - gf_var(Imm);
          grad_m = monocentral(var_slope_p, var_slope_m);
          var_slope_p = gf_var(Ipp) - gf_var(Ip);
          var_slope_m = gf_var(Ip) - gf_var(Im);
          grad_p = monocentral(var_slope_p, var_slope_m);
          right_face =
              (gf_var(Im) + gf_var(Ip)) / 2.0 + (grad_m - grad_p) / 6.0;

          // Now apply conditions in Eq. (1.11) of C&W (1984)
          qa = (right_face - gf_var(Im)) * (gf_var(Im) - left_face);
          qd = (right_face - left_face);
          qe = 6.0 * (gf_var(Im) - (left_face + right_face) / 2.0);
          if (qa <= 0.) {
            left_face = gf_var(Im);
            right_face = gf_var(Im);
          }
          if (qd * (qd - qe) < 0.0)
            left_face = 3.0 * gf_var(Im) - 2.0 * right_face;
          if (qd * (qd + qe) < 0.0)
            right_face = 3.0 * gf_var(Im) - 2.0 * left_face;

          // Keep right value of cell i-1 as i-1/2-eps
          CCTK_REAL var_m = right_face;

          return array<CCTK_REAL, 2>{var_m, var_p};
        }

        default:
          assert(0);
        }
      };

  const auto eigenvalues =
      [=] CCTK_DEVICE(CCTK_REAL alp_avg, CCTK_REAL beta_avg, CCTK_REAL u_avg,
                      array<CCTK_REAL, 2> vel, array<CCTK_REAL, 2> rho,
                      array<CCTK_REAL, 2> cs2, array<CCTK_REAL, 2> w_lor,
                      array<CCTK_REAL, 2> h) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // computing characteristics for the minus side
        array<CCTK_REAL, 3> a_m = {
            (cs2[0] * h[0] * rho[0]) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2[0]) * h[0] * rho[0] *
                    pow2(beta_avg - alp_avg * vel[0]) * pow2(w_lor[0]),

            2 * beta_avg * (cs2[0] * h[0] * rho[0]) -
                2 * (-1 + cs2[0]) * h[0] * rho[0] *
                    (beta_avg - alp_avg * vel[0]) * pow2(w_lor[0]),

            h[0] * rho[0] *
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
            (cs2[1] * h[1] * rho[1]) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2[1]) * h[1] * rho[1] *
                    pow2(beta_avg - alp_avg * vel[1]) * pow2(w_lor[1]),

            2 * beta_avg * (cs2[1] * h[1] * rho[1]) -
                2 * (-1 + cs2[1]) * h[1] * rho[1] *
                    (beta_avg - alp_avg * vel[1]) * pow2(w_lor[1]),

            h[1] * rho[1] *
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
    const array<CCTK_REAL, 2> rho_rc   = reconstruct(gf_rho,   p);
    const array<CCTK_REAL, 2> velx_rc  = reconstruct(gf_velx,  p);
    const array<CCTK_REAL, 2> vely_rc  = reconstruct(gf_vely,  p);
    const array<CCTK_REAL, 2> velz_rc  = reconstruct(gf_velz,  p);
    const array<CCTK_REAL, 2> eps_rc   = reconstruct(gf_eps,   p);
    const array<CCTK_REAL, 2> Bx_rc = reconstruct(gf_Bx, p);
    const array<CCTK_REAL, 2> By_rc = reconstruct(gf_By, p);
    const array<CCTK_REAL, 2> Bz_rc = reconstruct(gf_Bz, p);

    const array<array<CCTK_REAL, 2>, 3> vels_rc = {velx_rc, vely_rc, velz_rc};
    const array<CCTK_REAL, 2> vel_rc = vels_rc[dir];

    const array<array<CCTK_REAL, 2>, 3> Bs_rc = {Bx_rc, By_rc, Bz_rc};
    const array<CCTK_REAL, 2> B_rc = Bs_rc[dir];

    constexpr auto DI = PointDesc::DI;

    // TODO: to reconstruct w_lorentz*vel or 4-velocity u_i

    // Computing metric components
    CCTK_REAL alp_avg = 0;

    CCTK_REAL betax_avg = 0;
    CCTK_REAL betay_avg = 0;
    CCTK_REAL betaz_avg = 0;

    CCTK_REAL gxx_avg = 0;
    CCTK_REAL gxy_avg = 0;
    CCTK_REAL gxz_avg = 0;
    CCTK_REAL gyy_avg = 0;
    CCTK_REAL gyz_avg = 0;
    CCTK_REAL gzz_avg = 0;

    for (int dk = 0; dk < (dir == 2 ? 1 : 2); ++dk) {
      for (int dj = 0; dj < (dir == 1 ? 1 : 2); ++dj) {
        for (int di = 0; di < (dir == 0 ? 1 : 2); ++di) {
          alp_avg += gf_alp(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);

          betax_avg += gf_betax(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
          betay_avg += gf_betay(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
          betaz_avg += gf_betaz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);

          gxx_avg += gf_gxx(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
          gxy_avg += gf_gxy(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
          gxz_avg += gf_gxz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
          gyy_avg += gf_gyy(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
          gyz_avg += gf_gyz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
          gzz_avg += gf_gzz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
        }
      }
    }

    alp_avg /= 4;
    betax_avg /= 4;
    betay_avg /= 4;
    betaz_avg /= 4;
    gxx_avg /= 4;
    gxy_avg /= 4;
    gxz_avg /= 4;
    gyy_avg /= 4;
    gyz_avg /= 4;
    gzz_avg /= 4;

    const array<CCTK_REAL, 3> betas_avg = {betax_avg, betay_avg, betaz_avg};
    const CCTK_REAL beta_avg = betas_avg[dir];

    // TODO: Compute pressure based on user-specified EOS.
    // Currently, computing press for classical ideal gas from reconstructed
    // vars

    const array<CCTK_REAL, 2> press_rc = {eps_rc[0] * rho_rc[0] * (gamma - 1),
                                          eps_rc[1] * rho_rc[1] * (gamma - 1)};

    // Determinant of spatial metric
    const CCTK_REAL detg =
        calc_detg(gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg);
    const CCTK_REAL sqrt_detg = sqrt(detg);

    // Upper metric
    const array<CCTK_REAL, 6> ug_avg =
        calc_upperg(gxx_avg, gxy_avg, gxz_avg, gyy_avg, gyz_avg, gzz_avg, detg);

    // Array containing uxx, uyy, uzz components of the upper metric
    const array<CCTK_REAL, 3> ugs_avg = {ug_avg[0], ug_avg[2], ug_avg[5]};
    // Variable for either uxx, uyy or uzz depending on the direction
    const CCTK_REAL u_avg = ugs_avg[dir];

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


    // Computing the contravariant coordinate velocity v^i - beta^i/alpha using
    // the reconstructed variables
    const CCTK_REAL beta_over_alpha_avg = beta_avg / alp_avg;

    const array<CCTK_REAL, 2> velxshift_rc = {
        velx_rc[0] - beta_over_alpha_avg,
        velx_rc[1] - beta_over_alpha_avg};

    const array<CCTK_REAL, 2> velyshift_rc = {
        vely_rc[0] - beta_over_alpha_avg,
        vely_rc[1] - beta_over_alpha_avg};

    const array<CCTK_REAL, 2> velzshift_rc = {
        vely_rc[0] - beta_over_alpha_avg,
        vely_rc[1] - beta_over_alpha_avg};

    const array<array<CCTK_REAL, 2>, 3> velsshift_rc = {
        velxshift_rc, velyshift_rc, velzshift_rc};

    const array<CCTK_REAL, 2> velshift_rc = velsshift_rc[dir];


    // Computing the covariant coordinate velocity v^i - beta^i/alpha using the
    // reconstructed variables
    const array<CCTK_REAL, 2> velxshift_low_rc = {
        gxx_avg * velxshift_rc[0] + gxy_avg * velyshift_rc[0] + gxz_avg * velzshift_rc[0],
        gxx_avg * velxshift_rc[1] + gxy_avg * velyshift_rc[1] + gxz_avg * velzshift_rc[1]};

    const array<CCTK_REAL, 2> velyshift_low_rc = {
        gxy_avg * velxshift_rc[0] + gyy_avg * velyshift_rc[0] + gyz_avg * velzshift_rc[0],
        gxy_avg * velxshift_rc[1] + gyy_avg * velyshift_rc[1] + gyz_avg * velzshift_rc[1]};

    const array<CCTK_REAL, 2> velzshift_low_rc = {
        gxz_avg * velxshift_rc[0] + gyz_avg * velyshift_rc[0] + gzz_avg * velzshift_rc[0],
        gxz_avg * velxshift_rc[1] + gyz_avg * velyshift_rc[1] + gzz_avg * velzshift_rc[1]};



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
        Bx_rc[0]*Blowx_rc[0] + By_rc[0]*Blowy_rc[0] + Bz_rc[0]*Blowz_rc[0],
        Bx_rc[1]*Blowx_rc[1] + By_rc[1]*Blowy_rc[1] + Bz_rc[1]*Blowz_rc[1]};


    // Computing the magnetic field measured by the observer comoving with the
    // fluid using the reconstructed variables
    const array<CCTK_REAL, 2> alpha_b0_rc = {
        w_lorentz_rc[0]*(Bx_rc[0]*vlowx_rc[0] + By_rc[0]*vlowy_rc[0] + Bz_rc[0]*vlowz_rc[0]),
        w_lorentz_rc[1]*(Bx_rc[1]*vlowx_rc[1] + By_rc[0]*vlowy_rc[1] + Bz_rc[1]*vlowz_rc[1])};

    const array<CCTK_REAL, 2> blowx_rc = {
        Blowx_rc[0]/w_lorentz_rc[0] + alpha_b0_rc[0]*velxshift_low_rc[0],
        Blowx_rc[1]/w_lorentz_rc[1] + alpha_b0_rc[1]*velxshift_low_rc[1]};

    const array<CCTK_REAL, 2> blowy_rc = {
        Blowy_rc[0]/w_lorentz_rc[0] + alpha_b0_rc[0]*velyshift_low_rc[0],
        Blowy_rc[1]/w_lorentz_rc[1] + alpha_b0_rc[1]*velyshift_low_rc[1]};

    const array<CCTK_REAL, 2> blowz_rc = {
        Blowz_rc[0]/w_lorentz_rc[0] + alpha_b0_rc[0]*velzshift_low_rc[0],
        Blowz_rc[1]/w_lorentz_rc[1] + alpha_b0_rc[1]*velzshift_low_rc[1]};

    // FIXME: likely not needed by itself, but we'll see
    /*const array<CCTK_REAL, 2> b2_rc = {
        (B2_rc[0] + pow2(alpha_b0_rc[0]))/pow2(w_lorentz_rc[0]),
        (B2_rc[1] + pow2(alpha_b0_rc[1]))/pow2(w_lorentz_rc[1])};*/


    // Auxiliary variables to compute the conservative variables and their
    // fluxes
    const array<CCTK_REAL, 2> sqrt_detg_press_plus_pmag_rc = {
        sqrt_detg * (press_rc[0] + 0.5*(B2_rc[0] + pow2(alpha_b0_rc[0]))/pow2(w_lorentz_rc[0])),
        sqrt_detg * (press_rc[1] + 0.5*(B2_rc[1] + pow2(alpha_b0_rc[1]))/pow2(w_lorentz_rc[1]))};

    const array<CCTK_REAL, 2> sqrt_detg_B2_rc = {
        sqrt_detg * B2_rc[0],
        sqrt_detg * B2_rc[1]};

    const array<CCTK_REAL, 2> sqrt_detg_W2b2_rc = {
        sqrt_detg * (pow2(alpha_b0_rc[0]) + B2_rc[0]),
        sqrt_detg * (pow2(alpha_b0_rc[1]) + B2_rc[1])};

    const array<CCTK_REAL, 2> B_over_w_lorentz_rc = {
        B_rc[0] / w_lorentz_rc[0],
        B_rc[1] / w_lorentz_rc[1]};

    const array<CCTK_REAL, 2> alpha_b0_over_w_lorentz_rc = {
        alpha_b0_rc[0]/w_lorentz_rc[0],
        alpha_b0_rc[1]/w_lorentz_rc[1]};


    // Computing conservatives from primitives
    const array<CCTK_REAL, 2> dens_rc = {
        sqrt_detg * rho_rc[0] * w_lorentz_rc[0],
        sqrt_detg * rho_rc[1] * w_lorentz_rc[1]};

    const array<CCTK_REAL, 2> dens_h_W_rc = {
        dens_rc[0] * h_rc[0] * w_lorentz_rc[0],
        dens_rc[1] * h_rc[1] * w_lorentz_rc[1]};

    const array<CCTK_REAL, 2> dens_h_W_plus_sqrt_detg_W2b2_rc = {
        dens_h_W_rc[0] + sqrt_detg * (pow2(alpha_b0_rc[0]) + B2_rc[0]),
        dens_h_W_rc[1] + sqrt_detg * (pow2(alpha_b0_rc[1]) + B2_rc[1])};

    const array<CCTK_REAL, 2> momx_rc = {
        dens_h_W_plus_sqrt_detg_W2b2_rc[0] * vlowx_rc[0] - alpha_b0_rc[0]*blowx_rc[0],
        dens_h_W_plus_sqrt_detg_W2b2_rc[1] * vlowx_rc[1] - alpha_b0_rc[1]*blowx_rc[1]};

    const array<CCTK_REAL, 2> momy_rc = {
        dens_h_W_plus_sqrt_detg_W2b2_rc[0] * vlowy_rc[0] - alpha_b0_rc[0]*blowy_rc[0],
        dens_h_W_plus_sqrt_detg_W2b2_rc[1] * vlowy_rc[1] - alpha_b0_rc[1]*blowy_rc[1]};

    const array<CCTK_REAL, 2> momz_rc = {
        dens_h_W_plus_sqrt_detg_W2b2_rc[0] * vlowz_rc[0] - alpha_b0_rc[0]*blowz_rc[0],
        dens_h_W_plus_sqrt_detg_W2b2_rc[1] * vlowz_rc[1] - alpha_b0_rc[1]*blowz_rc[1]};

    // FIXME: B^2 = W^2·b^2 - (alpha·b^0)^2, is that true?
    const array<CCTK_REAL, 2> tau_rc = {
        dens_h_W_rc[0] - dens_rc[0] - sqrt_detg_press_plus_pmag_rc[0] + sqrt_detg_B2_rc[0],
        dens_h_W_rc[1] - dens_rc[1] - sqrt_detg_press_plus_pmag_rc[1] + sqrt_detg_B2_rc[1]};

    const array<CCTK_REAL, 2> Btildex_rc = {
        sqrt_detg * Bx_rc[0],
        sqrt_detg * Bx_rc[1]};

    const array<CCTK_REAL, 2> Btildey_rc = {
        sqrt_detg * By_rc[0],
        sqrt_detg * By_rc[1]};

    const array<CCTK_REAL, 2> Btildez_rc = {
        sqrt_detg * Bz_rc[0],
        sqrt_detg * Bz_rc[1]};


    // Computing fluxes of conserved variables
    const array<CCTK_REAL, 2> flux_dens = {dens_rc[0] * velshift_rc[0],
                                           dens_rc[1] * velshift_rc[1]};

    const array<CCTK_REAL, 2> flux_momx = {
        momx_rc[0] * velshift_rc[0] + (dir == 0) * sqrt_detg_press_plus_pmag_rc[0] - B_over_w_lorentz_rc[0] * blowx_rc[0],
        momx_rc[1] * velshift_rc[1] + (dir == 0) * sqrt_detg_press_plus_pmag_rc[1] - B_over_w_lorentz_rc[1] * blowx_rc[1]};

    const array<CCTK_REAL, 2> flux_momy = {
        momy_rc[0] * velshift_rc[0] + (dir == 1) * sqrt_detg_press_plus_pmag_rc[0] - B_over_w_lorentz_rc[0] * blowy_rc[0],
        momy_rc[1] * velshift_rc[1] + (dir == 1) * sqrt_detg_press_plus_pmag_rc[1] - B_over_w_lorentz_rc[1] * blowy_rc[1]};

    const array<CCTK_REAL, 2> flux_momz = {
        momz_rc[0] * velshift_rc[0] + (dir == 2) * sqrt_detg_press_plus_pmag_rc[0] - B_over_w_lorentz_rc[0] * blowz_rc[0],
        momz_rc[1] * velshift_rc[1] + (dir == 2) * sqrt_detg_press_plus_pmag_rc[1] - B_over_w_lorentz_rc[1] * blowz_rc[1]};

    const array<CCTK_REAL, 2> flux_tau = {
        tau_rc[0] * velshift_rc[0] + sqrt_detg_press_plus_pmag_rc[0] * vel_rc[0] - alpha_b0_over_w_lorentz_rc[0]*B_rc[0],
        tau_rc[1] * velshift_rc[1] + sqrt_detg_press_plus_pmag_rc[1] * vel_rc[1] - alpha_b0_over_w_lorentz_rc[1]*B_rc[1]};

    const array<CCTK_REAL, 2> flux_Btildex = {
        velshift_rc[0] * Btildex_rc[0] - velxshift_rc[0] * B_rc[0],
        velshift_rc[1] * Btildex_rc[1] - velxshift_rc[1] * B_rc[1]};

    const array<CCTK_REAL, 2> flux_Btildey = {
        velshift_rc[0] * Btildey_rc[0] - velyshift_rc[0] * B_rc[0],
        velshift_rc[1] * Btildey_rc[1] - velyshift_rc[1] * B_rc[1]};

    const array<CCTK_REAL, 2> flux_Btildez = {
        velshift_rc[0] * Btildez_rc[0] - velzshift_rc[0] * B_rc[0],
        velshift_rc[1] * Btildez_rc[1] - velzshift_rc[1] * B_rc[1]};


    array<array<CCTK_REAL, 4>, 2> lambda = eigenvalues(
        alp_avg, beta_avg, u_avg, vel_rc, rho_rc, cs2_rc, w_lorentz_rc, h_rc);


    gf_fluxdens(p.I) = calcflux(lambda, dens_rc, flux_dens);
    gf_fluxmomx(p.I) = calcflux(lambda, momx_rc, flux_momx);
    gf_fluxmomy(p.I) = calcflux(lambda, momy_rc, flux_momy);
    gf_fluxmomz(p.I) = calcflux(lambda, momz_rc, flux_momz);
    gf_fluxtau(p.I) = calcflux(lambda, tau_rc, flux_tau);
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

        const CCTK_REAL detg = calc_detg(gxx(p.I), gxy(p.I), gxz(p.I), gyy(p.I),
                                         gyz(p.I), gzz(p.I));
        const array<CCTK_REAL, 6> ug = calc_upperg(
            gxx(p.I), gxy(p.I), gxz(p.I), gyy(p.I), gyz(p.I), gzz(p.I), detg);

        const CCTK_REAL Axup =
            ug[0] * Ax_vert + ug[1] * Ay_vert + ug[2] * Az_vert;
        const CCTK_REAL Ayup =
            ug[1] * Ax_vert + ug[3] * Ay_vert + ug[4] * Az_vert;
        const CCTK_REAL Azup =
            ug[2] * Ax_vert + ug[4] * Ay_vert + ug[5] * Az_vert;

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
