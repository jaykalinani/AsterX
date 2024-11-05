#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include "c2p.hxx"
#include "c2p_1DPalenzuela.hxx"
#include "c2p_2DNoble.hxx"

#include "setup_eos.hxx"
#include "aster_utils.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace Con2PrimFactory;
using namespace AsterUtils;

enum class eos_3param { IdealGas, Hybrid, Tabulated };
enum class c2p_first_t { Noble, Palenzuela };
enum class c2p_second_t { Noble, Palenzuela };

template <typename EOSIDType, typename EOSType>
void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS, EOSIDType *eos_1p,
                             EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  c2p_first_t c2p_fir;
  c2p_second_t c2p_sec;

  if (CCTK_EQUALS(c2p_prime, "Noble")) {
    c2p_fir = c2p_first_t::Noble;
  } else if (CCTK_EQUALS(c2p_prime, "Palenzuela")) {
    c2p_fir = c2p_first_t::Palenzuela;
  } else {
    CCTK_ERROR("Unknown value for parameter \"c2p_prime\"");
  }

  if (CCTK_EQUALS(c2p_second, "Noble")) {
    c2p_sec = c2p_second_t::Noble;
  } else if (CCTK_EQUALS(c2p_second, "Palenzuela")) {
    c2p_sec = c2p_second_t::Palenzuela;
  } else {
    CCTK_ERROR("Unknown value for parameter \"c2p_second\"");
  }

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  // Loop over the interior of the grid
  cctk_grid.loop_int_device<
      1, 1, 1>(grid.nghostzones, [=] CCTK_DEVICE(
                                     const PointDesc
                                         &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    // Setting up atmosphere
    CCTK_REAL rho_atm = 0.0;   // dummy initialization
    CCTK_REAL press_atm = 0.0; // dummy initialization
    CCTK_REAL eps_atm = 0.0;   // dummy initialization
    CCTK_REAL local_eps_BH = eps_BH;
    CCTK_REAL radial_distance = sqrt(p.x * p.x + p.y * p.y + p.z * p.z);

    // Grading rho
    rho_atm = (radial_distance > r_atmo)
                  ? (rho_abs_min * pow((r_atmo / radial_distance), n_rho_atmo))
                  : rho_abs_min;
    const CCTK_REAL rho_atmo_cut = rho_atm * (1 + atmo_tol);

    // Grading pressure based on either cold or thermal EOS
    if (thermal_eos_atmo) {
      press_atm = (radial_distance > r_atmo)
                      ? (p_atmo * pow(r_atmo / radial_distance, n_press_atmo))
                      : p_atmo;
      eps_atm = eos_3p->eps_from_valid_rho_press_ye(rho_atm, press_atm, Ye_atmo);
    } else {
      const CCTK_REAL gm1 = eos_1p->gm1_from_valid_rho(rho_atm);
      eps_atm = eos_1p->sed_from_valid_gm1(gm1);
      eps_atm = std::min(std::max(eos_3p->rgeps.min, eps_atm), eos_3p->rgeps.max);
      press_atm = eos_3p->press_from_valid_rho_eps_ye(rho_atm, eps_atm, Ye_atmo);
    }
    atmosphere atmo(rho_atm, eps_atm, Ye_atmo, press_atm, rho_atmo_cut);

    // Construct Noble c2p object:
    c2p_2DNoble c2p_Noble(eos_3p, atmo, max_iter, c2p_tol, rho_strict, vw_lim,
                          B_lim, Ye_lenient);

    // Construct Palenzuela c2p object:
    c2p_1DPalenzuela c2p_Pal(eos_3p, atmo, max_iter, c2p_tol, rho_strict,
                             vw_lim, B_lim, Ye_lenient);

    /* Get covariant metric */
    const smat<CCTK_REAL, 3> glo(
        [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

    /* Calculate inverse of 3-metric */
    const CCTK_REAL spatial_detg = calc_det(glo);
    const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

    vec<CCTK_REAL, 3> v_up{saved_velx(p.I), saved_vely(p.I), saved_velz(p.I)};
    const vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);
    CCTK_REAL wlor = calc_wlorentz(v_low, v_up);

    vec<CCTK_REAL, 3> Bup{dBx(p.I) / sqrt_detg, dBy(p.I) / sqrt_detg,
                          dBz(p.I) / sqrt_detg};

    prim_vars pv;
    prim_vars pv_seeds{saved_rho(p.I), saved_eps(p.I), saved_Ye(p.I), press(p.I),
                       v_up,           wlor,           Bup};
    // Note that cv are densitized, i.e. they all include sqrt_detg
    cons_vars cv{dens(p.I),
                 {momx(p.I), momy(p.I), momz(p.I)},
                 tau(p.I),
                 DYe(p.I),
                 {dBx(p.I), dBy(p.I), dBz(p.I)}};

    if (dens(p.I) <= sqrt_detg * rho_atmo_cut) {
      cv.dBvec(0) = dBx(p.I); // densitized
      cv.dBvec(1) = dBy(p.I);
      cv.dBvec(2) = dBz(p.I);
      pv.Bvec = cv.dBvec / sqrt_detg;
      atmo.set(pv, cv, glo);
      atmo.set(pv_seeds);
    }

    // Modifying primitive seeds within BH interiors before C2Ps are called
    // NOTE: By default, alp_thresh=0 so the if condition below is never
    // triggered. One must be very careful when using this functionality and
    // must correctly set alp_thresh, rho_BH, eps_BH and vwlim_BH in the parfile

    if (alp(p.I) < alp_thresh) {
      if ((pv_seeds.rho > rho_BH) || (pv_seeds.eps > local_eps_BH)) {
        pv_seeds.rho = rho_BH; // typically set to 0.01% to 1% of rho_max of
                               // initial NS or disk
        pv_seeds.eps = local_eps_BH;
        pv_seeds.Ye = Ye_atmo;
        pv_seeds.press =
            eos_3p->press_from_valid_rho_eps_ye(rho_BH, local_eps_BH, Ye_atmo);
        // check on velocities
        CCTK_REAL wlim_BH = sqrt(1.0 + vwlim_BH * vwlim_BH);
        CCTK_REAL vlim_BH = vwlim_BH / wlim_BH;
        CCTK_REAL sol_v =
            sqrt((pv_seeds.w_lor * pv_seeds.w_lor - 1.0)) / pv_seeds.w_lor;
        if (sol_v > vlim_BH) {
          pv_seeds.vel *= vlim_BH / sol_v;
          pv_seeds.w_lor = wlim_BH;
        }
        cv.from_prim(pv_seeds, glo);
      }
    }

    // Construct error report object:
    c2p_report rep_first;
    c2p_report rep_second;

    // Calling the first C2P
    switch (c2p_fir) {
    case c2p_first_t::Noble: {
      c2p_Noble.solve(eos_3p, pv, pv_seeds, cv, glo, rep_first);
      break;
    }
    case c2p_first_t::Palenzuela: {
      c2p_Pal.solve(eos_3p, pv, pv_seeds, cv, glo, rep_first);
      break;
    }
    default:
      assert(0);
    }

    if (rep_first.failed()) {
      printf("First C2P failed :( \n");
      rep_first.debug_message();
      printf("Calling the back up C2P.. \n");
      // Calling the second C2P
      switch (c2p_sec) {
      case c2p_second_t::Noble: {
        c2p_Noble.solve(eos_3p, pv, pv_seeds, cv, glo, rep_second);
        break;
      }
      case c2p_second_t::Palenzuela: {
        c2p_Pal.solve(eos_3p, pv, pv_seeds, cv, glo, rep_second);
        break;
      }
      default:
        assert(0);
      }
    }

    if (rep_first.failed() && rep_second.failed()) {
      printf("Second C2P failed too :( :( \n");
      rep_second.debug_message();

      // Treatment for BH interiors after C2P failures
      // NOTE: By default, alp_thresh=0 so the if condition below is never
      // triggered. One must be very careful when using this functionality and
      // must correctly set alp_thresh, rho_BH, eps_BH and vwlim_BH in the
      // parfile
      if (alp(p.I) < alp_thresh) {
        if ((pv_seeds.rho > rho_BH) || (pv_seeds.eps > local_eps_BH)) {
          pv.rho = rho_BH; // typically set to 0.01% to 1% of rho_max of initial
                           // NS or disk
          pv.eps = local_eps_BH;
          pv.Ye = Ye_atmo;
          pv.press =
              eos_3p->press_from_valid_rho_eps_ye(rho_BH, local_eps_BH, Ye_atmo);
          // check on velocities
          CCTK_REAL wlim_BH = sqrt(1.0 + vwlim_BH * vwlim_BH);
          CCTK_REAL vlim_BH = vwlim_BH / wlim_BH;
          CCTK_REAL sol_v = sqrt((pv.w_lor * pv.w_lor - 1.0)) / pv.w_lor;
          if (sol_v > vlim_BH) {
            pv.vel *= vlim_BH / sol_v;
            pv.w_lor = wlim_BH;
          }
          cv.from_prim(pv, glo);
          rep_first.set_atmo = 0;
          rep_second.set_atmo = 0;
        }
      }
      con2prim_flag(p.I) = 0;
    }

    if (rep_first.set_atmo && rep_second.set_atmo) {
      if (debug_mode) {
        printf(
            "WARNING: \n"
            "C2Ps failed. Printing cons and saved prims before set to "
            "atmo: \n"
            "cctk_iteration = %i \n "
            "x, y, z = %26.16e, %26.16e, %26.16e \n "
            "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
            "momy = %26.16e \n momz = %26.16e \n DYe =  %26.16e \n dBx = %26.16e \n "
            "dBy = %26.16e \n dBz = %26.16e \n "
            "saved_rho = %26.16e \n saved_eps = %26.16e \n press= %26.16e \n "
            "saved_velx = %26.16e \n saved_vely = %26.16e \n saved_velz = "
            "%26.16e \n saved_Ye = %26.16e \n"
            "Bvecx = %26.16e \n Bvecy = %26.16e \n "
            "Bvecz = %26.16e \n "
            "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
            cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I), momx(p.I),
            momy(p.I), momz(p.I), DYe(p.I), dBx(p.I), dBy(p.I), dBz(p.I), pv.rho, pv.eps,
            pv.press, pv.vel(0), pv.vel(1), pv.vel(2), pv.Ye, pv.Bvec(0), pv.Bvec(1),
            pv.Bvec(2),
            // rho(p.I), eps(p.I), press(p.I), velx(p.I), vely(p.I),
            // velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I),
            Avec_x(p.I), Avec_y(p.I), Avec_z(p.I));
      }

      // set to atmo
      cv.dBvec(0) = dBx(p.I);
      cv.dBvec(1) = dBy(p.I);
      cv.dBvec(2) = dBz(p.I);
      pv.Bvec = cv.dBvec / sqrt_detg;
      atmo.set(pv, cv, glo);

      // assert(0);
    }

    /* set flag to success */
    con2prim_flag(p.I) = 1;

    // dummy vars
    CCTK_REAL Ex, Ey, Ez;

    // Write back pv
    pv.scatter(rho(p.I), eps(p.I), Ye(p.I), press(p.I), velx(p.I), vely(p.I),
               velz(p.I), wlor, Bvecx(p.I), Bvecy(p.I), Bvecz(p.I), Ex, Ey, Ez);

    zvec_x(p.I) = wlor * pv.vel(0);
    zvec_y(p.I) = wlor * pv.vel(1);
    zvec_z(p.I) = wlor * pv.vel(2);

    svec_x(p.I) =
        (pv.rho + pv.rho * pv.eps + pv.press) * wlor * wlor * pv.vel(0);
    svec_y(p.I) =
        (pv.rho + pv.rho * pv.eps + pv.press) * wlor * wlor * pv.vel(1);
    svec_z(p.I) =
        (pv.rho + pv.rho * pv.eps + pv.press) * wlor * wlor * pv.vel(2);

    // Write back cv
    cv.scatter(dens(p.I), momx(p.I), momy(p.I), momz(p.I), tau(p.I), DYe(p.I),
               dBx(p.I), dBy(p.I), dBz(p.I));

    // Update saved prims
    saved_rho(p.I) = rho(p.I);
    saved_velx(p.I) = velx(p.I);
    saved_vely(p.I) = vely(p.I);
    saved_velz(p.I) = velz(p.I);
    saved_eps(p.I) = eps(p.I);
    saved_Ye(p.I) = Ye(p.I);
  }); // Loop
}

extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  // defining EOS objects
  eos_3param eos_3p_type;

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_3p_type = eos_3param::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_3p_type = eos_3param::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_3p_type = eos_3param::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_3p_type) {
  case eos_3param::IdealGas: {
    // Get local eos objects
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_ig = global_eos_3p_ig;

    AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_1p_poly, eos_3p_ig);
    break;
  }
  case eos_3param::Hybrid: {
    // Get local eos objects
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_hyb = global_eos_3p_hyb;

    AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_1p_poly, eos_3p_hyb); 
    break;
  }
  case eos_3param::Tabulated: {
    // Get local eos objects
    auto eos_1p_poly = global_eos_1p_poly;
    auto eos_3p_tab3d = global_eos_3p_tab3d;

    AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_1p_poly, eos_3p_tab3d);
    break;
  }
  default:
    assert(0);
  }
}

extern "C" void AsterX_Con2Prim_Interpolate_Failed(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim_Interpolate_Failed;
  DECLARE_CCTK_PARAMETERS;

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const vec<GF3D2<CCTK_REAL>, 6> gf_prims{rho, velx, vely, velz, eps, press};
  const vec<GF3D2<CCTK_REAL>, 5> gf_cons{dens, momx, momy, momz, tau};

  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        if (con2prim_flag(p.I) == 0) {

          const vec<CCTK_REAL, 6> flag_nbs = get_neighbors(con2prim_flag, p);
          const vec<CCTK_REAL, 6> rho_nbs = get_neighbors(rho, p);
          const vec<CCTK_REAL, 6> velx_nbs = get_neighbors(velx, p);
          const vec<CCTK_REAL, 6> vely_nbs = get_neighbors(vely, p);
          const vec<CCTK_REAL, 6> velz_nbs = get_neighbors(velz, p);
          const vec<CCTK_REAL, 6> eps_nbs = get_neighbors(eps, p);
          const vec<CCTK_REAL, 6> saved_rho_nbs = get_neighbors(saved_rho, p);
          const vec<CCTK_REAL, 6> saved_velx_nbs = get_neighbors(saved_velx, p);
          const vec<CCTK_REAL, 6> saved_vely_nbs = get_neighbors(saved_vely, p);
          const vec<CCTK_REAL, 6> saved_velz_nbs = get_neighbors(saved_velz, p);
          const vec<CCTK_REAL, 6> saved_eps_nbs = get_neighbors(saved_eps, p);

          CCTK_REAL sum_nbs =
              sum<6>([&](int i) ARITH_INLINE { return flag_nbs(i); });
          assert(sum_nbs > 0);
          rho(p.I) = calc_avg_neighbors(flag_nbs, rho_nbs, saved_rho_nbs);
          velx(p.I) = calc_avg_neighbors(flag_nbs, velx_nbs, saved_velx_nbs);
          vely(p.I) = calc_avg_neighbors(flag_nbs, vely_nbs, saved_vely_nbs);
          velz(p.I) = calc_avg_neighbors(flag_nbs, velz_nbs, saved_velz_nbs);
          eps(p.I) = calc_avg_neighbors(flag_nbs, eps_nbs, saved_eps_nbs);
          press(p.I) = (gl_gamma - 1) * eps(p.I) * rho(p.I);

          /* reset flag */
          con2prim_flag(p.I) = 1;

          // set to atmos
          /*
          if (rho(p.I) <= rho_abs_min * (1 + atmo_tol))
          {
            const smat<CCTK_REAL, 3> g3_avg([&](int i, int j) ARITH_INLINE
                                            { return calc_avg_v2c(gf_g(i, j),
          p); }); const CCTK_REAL sqrtg = sqrt(calc_det(g3_avg)); const
          vec<CCTK_REAL, 3> Bup{Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)}; const
          vec<CCTK_REAL, 3> Blow = calc_contraction(g3_avg, Bup); const
          CCTK_REAL Bsq = calc_contraction(Bup, Blow);

            set_to_atmosphere(rho_abs_min, poly_K, gamma, sqrtg, Bsq, gf_prims,
                              gf_cons, p);
          };
          */
          saved_rho(p.I) = rho(p.I);
          saved_velx(p.I) = velx(p.I);
          saved_vely(p.I) = vely(p.I);
          saved_velz(p.I) = velz(p.I);
          saved_eps(p.I) = eps(p.I);
        }
      });
}

} // namespace AsterX
