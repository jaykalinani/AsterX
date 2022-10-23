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
                      vec<CCTK_REAL, 2> vel, vec<CCTK_REAL, 2> rho,
                      vec<CCTK_REAL, 2> cs2, vec<CCTK_REAL, 2> w_lor,
                      vec<CCTK_REAL, 2> h,
                      vec<CCTK_REAL, 2> bsq) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // computing characteristics for the minus side
        // See Eq. (28) of Giacomazzo & Rezzolla (2007) with b^i=0
        array<CCTK_REAL, 3> a_m = {
            (bsq(0) + cs2(0) * h(0) * rho(0)) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2(0)) * h(0) * rho(0) *
                    pow2(beta_avg - alp_avg * vel(0)) * pow2(w_lor(0)),

            2 * beta_avg * (bsq(0) + cs2(0) * h(0) * rho(0)) -
                2 * (-1 + cs2(0)) * h(0) * rho(0) *
                    (beta_avg - alp_avg * vel(0)) * pow2(w_lor(0)),

            bsq(0) + h(0) * rho(0) *
                         (cs2(0) + pow2(w_lor(0)) - cs2(0) * pow2(w_lor(0)))};

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
            (bsq(1) + cs2(1) * h(1) * rho(1)) *
                    (pow2(beta_avg) - pow2(alp_avg) * u_avg) -
                (-1 + cs2(1)) * h(1) * rho(1) *
                    pow2(beta_avg - alp_avg * vel(1)) * pow2(w_lor(1)),

            2 * beta_avg * (bsq(1) + cs2(1) * h(1) * rho(1)) -
                2 * (-1 + cs2(1)) * h(1) * rho(1) *
                    (beta_avg - alp_avg * vel(1)) * pow2(w_lor(1)),

            bsq(1) + h(1) * rho(1) *
                         (cs2(1) + pow2(w_lor(1)) - cs2(1) * pow2(w_lor(1)))};

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
      [=] CCTK_DEVICE(array<array<CCTK_REAL, 4>, 2> lam, vec<CCTK_REAL, 2> var,
                      vec<CCTK_REAL, 2> flux) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        const CCTK_REAL charmax =
            max({0.0, fabs(lam[0][0]), fabs(lam[0][1]), fabs(lam[0][2]),
                 fabs(lam[0][3]), fabs(lam[1][0]), fabs(lam[1][1]),
                 fabs(lam[1][2]), fabs(lam[1][3])});

        CCTK_REAL llf =
            0.5 * ((flux(0) + flux(1)) - charmax * (var(1) - var(0)));
        // return dA * llf;
        return llf;
      };

  grid.loop_int_device<face_centred[0], face_centred[1], face_centred[2]>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        /* vec of grid functions */
        const vec<GF3D2<const CCTK_REAL>, 3> gf_vels{velx, vely, velz};
        const vec<GF3D2<const CCTK_REAL>, 3> gf_Bvecs{Bvecx, Bvecy, Bvecz};
        const vec<GF3D2<const CCTK_REAL>, 3> gf_beta{betax, betay, betaz};
        const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz,
                                                   gyy, gyz, gzz};

        // Reconstruct primitives from the cells on left (indice 0) and right
        // (indice 1) side of this face rc = reconstructed variables or computed
        // from reconstructed variables

        const vec<CCTK_REAL, 2> rho_rc{reconstruct_pt(rho, p)};
        const vec<vec<CCTK_REAL, 2>, 3> vels_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_pt(gf_vels(i), p)};
        });
        const vec<CCTK_REAL, 2> eps_rc{reconstruct_pt(eps, p)};
        const vec<vec<CCTK_REAL, 2>, 3> Bs_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>{reconstruct_pt(gf_Bvecs(i), p)};
        });

        // Computing metric components
        const CCTK_REAL alp_avg = calc_avg_v2f(alp, p, dir);
        const vec<CCTK_REAL, 3> betas_avg([&](int i) ARITH_INLINE {
          return calc_avg_v2f(gf_beta(i), p, dir);
        });
        const smat<CCTK_REAL, 3> g_avg([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2f(gf_g(i, j), p, dir);
        });
        // Determinant of spatial metric
        const CCTK_REAL detg_avg = calc_det(g_avg);
        const CCTK_REAL sqrt_detg = sqrt(detg_avg);
        // Upper metric
        const smat<CCTK_REAL, 3> ug_avg = calc_inv(g_avg, detg_avg);
        // Variable for either uxx, uyy or uzz depending on the direction
        const CCTK_REAL u_avg = ug_avg(dir, dir);

        // TODO: Compute pressure based on user-specified EOS.
        // Currently, computing press for classical ideal gas from reconstructed
        // vars
        const vec<CCTK_REAL, 2> press_rc([&](int f) ARITH_INLINE {
          return eps_rc(f) * rho_rc(f) * (gamma - 1);
        });

        // v_j
        const vec<vec<CCTK_REAL, 2>, 3> vlows_rc =
            calc_contraction(g_avg, vels_rc);
        // Computing the contravariant coordinate velocity
        // vtilde^i = alpha*v^i - beta^i using the reconstructed variables
        // const vec<vec<CCTK_REAL, 3>, 2> vtilde_rc = alp_avg * vel_c -
        // beta_avg;
        const vec<vec<CCTK_REAL, 2>, 3> vtildes_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return alp_avg * vels_rc(i)(f) - betas_avg(i);
          });
        });

        // Computing w_lorentz using reconstructed variables
        const vec<CCTK_REAL, 2> w_lorentz_rc([&](int f) ARITH_INLINE {
          return 1.0 / sqrt(1.0 - calc_contraction(vlows_rc, vels_rc)(f));
        });

        // Computing cs2 for ideal gas EOS using reconstructed variables
        const vec<CCTK_REAL, 2> cs2_rc([&](int f) ARITH_INLINE {
          return (gamma - 1.0) * eps_rc(f) / (eps_rc(f) + 1.0 / gamma);
        });

        // Computing enthalpy h for ideal gas EOS using reconstructed variables
        const vec<CCTK_REAL, 2> h_rc([&](int f) ARITH_INLINE {
          return 1.0 + eps_rc(f) + press_rc(f) / rho_rc(f);
        });

        // Computing the covariant magnetic field measured by the Eulerian
        // observer using the reconstructed variables
        const vec<vec<CCTK_REAL, 2>, 3> Blows_rc =
            calc_contraction(g_avg, Bs_rc);
        const vec<CCTK_REAL, 2> B2_rc = calc_contraction(Bs_rc, Blows_rc);

        // Computing the magnetic field measured by the observer comoving with
        // the fluid using the reconstructed variables
        const vec<CCTK_REAL, 2> alpha_b0_rc([&](int f) ARITH_INLINE {
          return w_lorentz_rc(f) * calc_contraction(Bs_rc, vlows_rc)(f);
        });

        const vec<vec<CCTK_REAL, 2>, 3> blows_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return Blows_rc(i)(f) / w_lorentz_rc(f) +
                   alpha_b0_rc(f) * vlows_rc(i)(f);
          });
        });

        const vec<CCTK_REAL, 2> bsq_rc([&](int f) ARITH_INLINE {
          return (B2_rc(f) + pow2(alpha_b0_rc(f))) / pow2(w_lorentz_rc(f));
        });

        /* Componets correspond to the dir we are considering */
        const CCTK_REAL beta_avg = betas_avg(dir);
        const vec<CCTK_REAL, 2> vel_rc{vels_rc(dir)};
        const vec<CCTK_REAL, 2> B_rc{Bs_rc(dir)};
        const vec<CCTK_REAL, 2> vtilde_rc{vtildes_rc(dir)};

        // Auxiliary variables to compute the conservative variables and their
        // fluxes
        const vec<CCTK_REAL, 2> sqrt_detg_press_plus_pmag_rc =
            sqrt_detg * (press_rc + 0.5 * bsq_rc);
        const vec<CCTK_REAL, 2> alp_sqrt_detg_press_plus_pmag_rc =
            alp_avg * sqrt_detg_press_plus_pmag_rc;

        const vec<CCTK_REAL, 2> alp_sqrt_detg_B_over_w_lorentz_rc(
            [&](int f) ARITH_INLINE {
              return alp_avg * sqrt_detg * B_rc(f) / w_lorentz_rc(f);
            });

        // Computing conservatives from primitives
        const vec<CCTK_REAL, 2> dens_rc([&](int f) ARITH_INLINE {
          return sqrt_detg * rho_rc(f) * w_lorentz_rc(f);
        });

        const vec<CCTK_REAL, 2> dens_h_W_rc([&](int f) ARITH_INLINE {
          return dens_rc(f) * h_rc(f) * w_lorentz_rc(f);
        });

        // sqrt(g)*( rho*h*W^2 + b^2*W^2 ) = sqrt(g)( rho*h*W^2 + (alp^2*b0^2) +
        // B^2 )
        const vec<CCTK_REAL, 2> dens_h_W_plus_sqrt_detg_W2b2_rc(
            [&](int f) ARITH_INLINE {
              return dens_h_W_rc(f) +
                     sqrt_detg * (pow2(alpha_b0_rc(f)) + B2_rc(f));
            });

        const vec<vec<CCTK_REAL, 2>, 3> moms_rc([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return dens_h_W_plus_sqrt_detg_W2b2_rc(f) * vlows_rc(i)(f) -
                   sqrt_detg * alpha_b0_rc(f) * blows_rc(i)(f);
          });
        });

        const vec<CCTK_REAL, 2> tau_rc = dens_h_W_rc - dens_rc -
                                         sqrt_detg_press_plus_pmag_rc +
                                         sqrt_detg * B2_rc;

        const vec<vec<CCTK_REAL, 2>, 3> Btildes_rc(
            [&](int i) ARITH_INLINE { return sqrt_detg * Bs_rc(i); });

        // Computing fluxes of conserved variables
        const vec<CCTK_REAL, 2> flux_dens(
            [&](int f) ARITH_INLINE { return dens_rc(f) * vtilde_rc(f); });

        const vec<CCTK_REAL, 3> vdirs{(dir == 0), (dir == 1), (dir == 2)};
        const vec<vec<CCTK_REAL, 2>, 3> flux_moms([&](int i) ARITH_INLINE {
          return vec<CCTK_REAL, 2>([&](int f) ARITH_INLINE {
            return moms_rc(i)(f) * vtilde_rc(f) +
                   vdirs(i) * alp_sqrt_detg_press_plus_pmag_rc(f) -
                   alp_sqrt_detg_B_over_w_lorentz_rc(f) * blows_rc(i)(f);
          });
        });

        const vec<CCTK_REAL, 2> flux_tau([&](int f) ARITH_INLINE {
          return tau_rc(f) * vtilde_rc(f) +
                 alp_sqrt_detg_press_plus_pmag_rc(f) * vel_rc(f) -
                 alpha_b0_rc(f) * alp_sqrt_detg_B_over_w_lorentz_rc(f);
        });

        // (0, Ez, -Ey), (-Ez, 0, Ex), (Ey, -Ex, 0)
        const vec<vec<CCTK_REAL, 2>, 3> Es_rc =
            calc_cross_product(Btildes_rc, vtildes_rc);
        const vec<vec<CCTK_REAL, 2>, 3> flux_Btildes =
            calc_cross_product(vdirs, Es_rc);

        array<array<CCTK_REAL, 4>, 2> lambda =
            eigenvalues(alp_avg, beta_avg, u_avg, vel_rc, rho_rc, cs2_rc,
                        w_lorentz_rc, h_rc, bsq_rc);

        fluxdenss[dir](p.I) = calcflux(lambda, dens_rc, flux_dens);
        fluxmomxs[dir](p.I) = calcflux(lambda, moms_rc(0), flux_moms(0));
        fluxmomys[dir](p.I) = calcflux(lambda, moms_rc(1), flux_moms(1));
        fluxmomzs[dir](p.I) = calcflux(lambda, moms_rc(2), flux_moms(2));
        fluxtaus[dir](p.I) = calcflux(lambda, tau_rc, flux_tau);
        fluxBxs[dir](p.I) =
            (dir != 0) * calcflux(lambda, Btildes_rc(0), flux_Btildes(0));
        fluxBys[dir](p.I) =
            (dir != 1) * calcflux(lambda, Btildes_rc(1), flux_Btildes(1));
        fluxBzs[dir](p.I) =
            (dir != 2) * calcflux(lambda, Btildes_rc(2), flux_Btildes(2));
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
            ug(0, 0) * Ax_vert + ug(0, 1) * Ay_vert + ug(0, 2) * Az_vert;
        const CCTK_REAL Ayup =
            ug(0, 1) * Ax_vert + ug(1, 1) * Ay_vert + ug(1, 2) * Az_vert;
        const CCTK_REAL Azup =
            ug(0, 2) * Ax_vert + ug(1, 2) * Ay_vert + ug(2, 2) * Az_vert;

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
