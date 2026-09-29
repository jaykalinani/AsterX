#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include "c2p.hxx"
#include "c2p_1DEntropy.hxx"
#include "c2p_1DPalenzuela.hxx"
#include "c2p_1DRePrimAnd.hxx"
#include "c2p_2DNoble.hxx"

#include "aster_utils.hxx"
#include "setup_eos.hxx"
#include "repair_diagnostics.hxx"

namespace AsterX {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace Con2PrimFactory;
using namespace AsterUtils;

enum class eos_3param { IdealGas, Hybrid, Tabulated };
enum class c2p_first_t { None, Noble, RePrimAnd, Palenzuela, Entropy };
enum class c2p_second_t { None, Noble, RePrimAnd, Palenzuela, Entropy };

enum C2PFlag : CCTK_INT {
  C2P_INIT = 0,    // initial value
  C2P_PRIME = 1,   // first solver succeeded
  C2P_SECOND = 2,  // second solver succeeded
  C2P_ENTROPY = 3, // 1‑D Entropy (kappa) solver succeeded
  C2P_ATMO = 4,    // conservative density lies below the atmosphere cutoff
  C2P_AVG = 5,     // primitives obtained by neighbour‑averaging
  C2P_FAIL = 6     // when C2P fails
};

template <typename EOSIDType, typename EOSType>
void AsterX_Con2Prim_typeEoS(CCTK_ARGUMENTS, EOSIDType *eos_1p,
                             EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  repair_diagnostics diagnostics(repair_every > 0 &&
      cctk_iteration % repair_every == 0);
  auto *counts = diagnostics.data();

  c2p_first_t c2p_fir;
  c2p_second_t c2p_sec;

  if (CCTK_EQUALS(c2p_prime, "Noble")) {
    c2p_fir = c2p_first_t::Noble;
  } else if (CCTK_EQUALS(c2p_prime, "RePrimAnd")) {
    c2p_fir = c2p_first_t::RePrimAnd;
  } else if (CCTK_EQUALS(c2p_prime, "Palenzuela")) {
    c2p_fir = c2p_first_t::Palenzuela;
  } else if (CCTK_EQUALS(c2p_prime, "Entropy")) {
    c2p_fir = c2p_first_t::Entropy;
  } else if (CCTK_EQUALS(c2p_prime, "None")) {
    c2p_fir = c2p_first_t::None;
  } else {
    CCTK_ERROR("Unknown value for parameter \"c2p_prime\"");
  }

  if (CCTK_EQUALS(c2p_second, "Noble")) {
    c2p_sec = c2p_second_t::Noble;
  } else if (CCTK_EQUALS(c2p_second, "RePrimAnd")) {
    c2p_sec = c2p_second_t::RePrimAnd;
  } else if (CCTK_EQUALS(c2p_second, "Palenzuela")) {
    c2p_sec = c2p_second_t::Palenzuela;
  } else if (CCTK_EQUALS(c2p_second, "Entropy")) {
    c2p_sec = c2p_second_t::Entropy;
  } else if (CCTK_EQUALS(c2p_second, "None")) {
    c2p_sec = c2p_second_t::None;
  } else {
    CCTK_ERROR("Unknown value for parameter \"c2p_second\"");
  }

  const vec<GF3D2<const CCTK_REAL>, dim> gf_beta{betax, betay, betaz};

  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};

  const auto c2p_impl = [=] CCTK_DEVICE(
                            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
    // Note that HydroBaseX gfs are NaN when entering this loop due
    // explicit dependence on conservatives from
    // AsterX -> dependents tag

    const CCTK_REAL radial_distance =
        sqrt(p.x * p.x + p.y * p.y + p.z * p.z);
    const auto atmo = make_atmosphere(
        eos_1p, eos_3p, radial_distance, rho_abs_min, p_atmo, t_atmo,
        Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo, n_temp_atmo, atmo_tol,
        thermal_eos_atmo, use_press_atmo);

    // ----- Construct C2P objects -----

    // Construct Noble c2p object:
    c2p_2DNoble c2p_Noble(eos_3p, atmo, max_iter, c2p_tol, alp_thresh, vw_lim,
                          B_lim, rho_BH, eps_BH, vwlim_BH, sigma_max,
                          inv_beta_max, Ye_lenient, use_z, use_temperature,
                          use_press_atmo, soft_root_convergence,
                          soft_root_width_factor);

    // Construct RePrimAnd c2p object:
    c2p_1DRePrimAnd c2p_RPA(eos_3p, atmo, max_iter, c2p_tol, alp_thresh, vw_lim,
                            B_lim, rho_BH, eps_BH, vwlim_BH, sigma_max,
                            inv_beta_max, Ye_lenient, use_z, use_temperature,
                            use_press_atmo, soft_root_convergence,
                            soft_root_width_factor);

    // Construct Palenzuela c2p object:
    c2p_1DPalenzuela c2p_Pal(eos_3p, atmo, max_iter, c2p_tol, alp_thresh,
                             vw_lim, B_lim, rho_BH, eps_BH, vwlim_BH, sigma_max,
                             inv_beta_max, Ye_lenient, use_z, use_temperature,
                             use_press_atmo, soft_root_convergence,
                             soft_root_width_factor);

    // Construct Entropy c2p object:
    c2p_1DEntropy c2p_Ent(eos_3p, atmo, max_iter, c2p_tol, alp_thresh, vw_lim,
                          B_lim, rho_BH, eps_BH, vwlim_BH, sigma_max,
                          inv_beta_max, Ye_lenient, use_z, use_temperature,
                          use_press_atmo, soft_root_convergence,
                          soft_root_width_factor);

    // ----------

    /* Get lapse */
    const CCTK_REAL alp_avg = calc_avg_v2c(alp, p);

    /* Get shift */
    const vec<CCTK_REAL, 3> beta_avg(
        [&](int i) ARITH_INLINE { return calc_avg_v2c(gf_beta(i), p); });

    /* Get covariant metric */
    const smat<CCTK_REAL, 3> glo(
        [&](int i, int j) ARITH_INLINE { return calc_avg_v2c(gf_g(i, j), p); });

    /* Get mask */
    CCTK_REAL mask_local = use_mask ? calc_avg_v2c(aster_mask_vc, p) : 1.0;

    /* Calculate inverse of 3-metric */
    const CCTK_REAL spatial_detg = calc_det(glo);
    const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

    vec<CCTK_REAL, 3> v_up{saved_velx(p.I), saved_vely(p.I), saved_velz(p.I)};
    vec<CCTK_REAL, 3> v_low = calc_contraction(glo, v_up);
    CCTK_REAL zsq{0.0};

    // TODO: Debug code to capture v>1 early,
    // remove soon
    const CCTK_REAL vsq = calc_contraction(v_low, v_up);
    if (vsq >= 1.0) {
      CCTK_REAL wlim = sqrt(1.0 + vw_lim * vw_lim);
      CCTK_REAL vlim = vw_lim / wlim;
      v_up *= vlim / sqrt(vsq);
      v_low *= vlim / sqrt(vsq);
      zsq = vw_lim * vw_lim;
    } else {
      zsq = vsq / (1.0 - vsq);
    }

    // CCTK_REAL wlor = calc_wlorentz(v_low, v_up);
    CCTK_REAL wlor = sqrt(1.0 + zsq);

    // Note that cv are densitized, i.e. they all include sqrt_detg
    cons_vars cv{dens(p.I), {momx(p.I), momy(p.I), momz(p.I)},
                 tau(p.I),  DYe(p.I),
                 DEnt(p.I), {dBx(p.I), dBy(p.I), dBz(p.I)}};

    // Undensitized magnetic fields
    const vec<CCTK_REAL, 3> Bup{cv.dBvec(0) / sqrt_detg,
                                cv.dBvec(1) / sqrt_detg,
                                cv.dBvec(2) / sqrt_detg};

    prim_vars pv;
    prim_vars pv_seeds{saved_rho(p.I),
                       saved_eps(p.I),
                       saved_Ye(p.I),
                       eos_3p->press_from_rho_eps_ye(
                           saved_rho(p.I), saved_eps(p.I), saved_Ye(p.I)),
                       temperature(p.I),
                       eos_3p->kappa_from_rho_eps_ye(
                           saved_rho(p.I), saved_eps(p.I), saved_Ye(p.I)),
                       v_up,
                       wlor,
                       Bup};

    /* set flag to success */
    bool c2p_flag_local = true;
    CCTK_INT c2p_flag_code = C2P_INIT;
    bool call_c2p = true;

    // Check if point is below atmosphere, and if atmosphere obeys magnetization
    // limits (RPA only). Magnetization limits are currently only applied for RPA C2P, 
    // while they are not obeyed in the other cases in the atmopshere -> TODO
    const CCTK_REAL b2_atm = calc_norm(Bup, glo);
    const bool set_atmo = (cv.dens <= sqrt_detg * atmo.rho_cut) &&
                          (c2p_off_floor_strict ||
                           ((b2_atm / atmo.rho_atmo <= sigma_max) &&
                            (b2_atm / (2 * atmo.press_atmo) <= inv_beta_max)));
    if (set_atmo) {
      pv.Bvec = Bup;
      atmo.set(pv, cv, glo);
      atmo.set(pv_seeds);
      c2p_flag_code = C2P_ATMO;
      call_c2p = false;
    }

    // Modifying primitive seeds within BH interiors before C2Ps are called
    // NOTE: By default, alp_thresh=0 so the if condition below is never
    // triggered. One must be very careful when using this functionality and
    // must correctly set alp_thresh, rho_BH, eps_BH and vwlim_BH in the parfile

    if (alp_avg < alp_thresh) {
      mask_local = 0.0;
    }
    aster_mask_cc(p.I) = mask_local;

    if (excise) {

      if (mask_local != 1.0) {
        c2p_Noble.bh_interior<EOSType, false>(eos_3p, pv_seeds, cv, glo);
        pv = pv_seeds;
        call_c2p = false;
      }
    }

    // Construct error report object:
    c2p_report rep_first;
    c2p_report rep_second;
    c2p_report rep_ent;

    // ----- ----- C2P ----- -----

    if (call_c2p) {

      // Limit conservatives before calling C2P
      const CCTK_REAL tau_before = cv.tau;
      c2p_Noble.cons_floors_and_ceilings(eos_3p, cv, glo, tauFluid_atmo);
      count_repair(counts, tau_repair, cv.tau != tau_before);

      // Calling the first C2P
      c2p_flag_code = C2P_PRIME;
      switch (c2p_fir) {
      case c2p_first_t::Noble: {
        c2p_Noble.solve(eos_3p, pv, pv_seeds, cv, alp_avg, beta_avg, glo,
                        rep_first, use_entropy_fix);
        break;
      }
      case c2p_first_t::RePrimAnd: {
        c2p_RPA.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_first,
                      use_entropy_fix);
        break;
      }
      case c2p_first_t::Palenzuela: {
        c2p_Pal.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_first,
                      use_entropy_fix);
        break;
      }
      case c2p_first_t::Entropy: {
        c2p_Ent.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_first);
        break;
      }
      case c2p_first_t::None: {
        // solve not called, pv remains unwritten
        break;
      }
      default:
        assert(0);
      }

      if (rep_first.failed()) {
        count_repair(counts, primary_failure,
                     c2p_fir != c2p_first_t::None);
        count_repair(counts, backup_call,
                     c2p_sec != c2p_second_t::None);
        c2p_flag_code = C2P_SECOND;

        if (debug_mode) {
          printf("First C2P failed :( \n");
          rep_first.debug_message();
          printf("Calling the back up C2P.. \n");
        }

        // Calling the second C2P
        switch (c2p_sec) {
        case c2p_second_t::Noble: {
          c2p_Noble.solve(eos_3p, pv, pv_seeds, cv, alp_avg, beta_avg, glo,
                          rep_second, use_entropy_fix);
          break;
        }
        case c2p_second_t::RePrimAnd: {
          c2p_RPA.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_second,
                        use_entropy_fix);
          break;
        }
        case c2p_second_t::Palenzuela: {
          c2p_Pal.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_second,
                        use_entropy_fix);
          break;
        }
        case c2p_second_t::Entropy: {
          c2p_Ent.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_second);
          break;
        }
        case c2p_second_t::None: {
          // solve not called, pv remains unwritten
          break;
        }
        default:
          assert(0);
        }
      }

      if (rep_first.failed() && rep_second.failed()) {

        if (use_entropy_fix) {

          c2p_flag_code = C2P_ENTROPY;
          c2p_Ent.solve(eos_3p, pv, cv, alp_avg, beta_avg, glo, rep_ent);

          if (rep_ent.failed()) {

            c2p_flag_local = false;
            c2p_flag_code = C2P_FAIL;

            if (debug_mode) {
              printf("Entropy C2P failed. Setting point to atmosphere.\n");
              rep_ent.debug_message();
              printf("WARNING: \n"
                     "C2Ps failed. Printing cons and saved prims before set to "
                     "atmo: \n"
                     "cctk_iteration = %i \n "
                     "x, y, z = %26.16e, %26.16e, %26.16e \n "
                     "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
                     "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
                     "dBy = %26.16e \n dBz = %26.16e \n "
                     "saved_rho = %26.16e \n saved_eps = %26.16e \n press= "
                     "%26.16e "
                     "\n "
                     "saved_velx = %26.16e \n saved_vely = %26.16e \n "
                     "saved_velz = "
                     "%26.16e \n "
                     "Bvecx = %26.16e \n Bvecy = %26.16e \n "
                     "Bvecz = %26.16e \n "
                     "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e "
                     "\n ",
                     cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I),
                     momx(p.I), momy(p.I), momz(p.I), dBx(p.I), dBy(p.I),
                     dBz(p.I), pv.rho, pv.eps, pv.press, pv.vel(0), pv.vel(1),
                     pv.vel(2), pv.Bvec(0), pv.Bvec(1), pv.Bvec(2),
                     // rho(p.I), eps(p.I), press(p.I), velx(p.I), vely(p.I),
                     // velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I),
                     Avec_x(p.I), Avec_y(p.I), Avec_z(p.I));
            }

            if (mask_local != 1.0) {
              // Failure inside mask
              c2p_Noble.bh_interior<EOSType, false>(eos_3p, pv_seeds, cv, glo);
              pv = pv_seeds;
            } else {
              // Failure outside, set to atmo
              cv.dBvec(0) = sqrt_detg * Bup(0);
              cv.dBvec(1) = sqrt_detg * Bup(1);
              cv.dBvec(2) = sqrt_detg * Bup(2);
              pv.Bvec = Bup;
              atmo.set(pv, cv, glo);
            }
          }

        } else {

          c2p_flag_local = false;
          c2p_flag_code = C2P_FAIL;

          if (debug_mode) {
            printf("Second C2P failed too :( :( \n");
            rep_second.debug_message();
            printf(
                "WARNING: \n"
                "C2Ps failed. Printing cons and saved prims before set to "
                "atmo: \n"
                "cctk_iteration = %i \n "
                "x, y, z = %26.16e, %26.16e, %26.16e \n "
                "dens = %26.16e \n tau = %26.16e \n momx = %26.16e \n "
                "momy = %26.16e \n momz = %26.16e \n dBx = %26.16e \n "
                "dBy = %26.16e \n dBz = %26.16e \n "
                "saved_rho = %26.16e \n saved_eps = %26.16e \n press= %26.16e "
                "\n "
                "saved_velx = %26.16e \n saved_vely = %26.16e \n saved_velz = "
                "%26.16e \n "
                "Bvecx = %26.16e \n Bvecy = %26.16e \n "
                "Bvecz = %26.16e \n "
                "Avec_x = %26.16e \n Avec_y = %26.16e \n Avec_z = %26.16e \n ",
                cctk_iteration, p.x, p.y, p.z, dens(p.I), tau(p.I), momx(p.I),
                momy(p.I), momz(p.I), dBx(p.I), dBy(p.I), dBz(p.I), pv.rho,
                pv.eps, pv.press, pv.vel(0), pv.vel(1), pv.vel(2), pv.Bvec(0),
                pv.Bvec(1), pv.Bvec(2),
                // rho(p.I), eps(p.I), press(p.I), velx(p.I), vely(p.I),
                // velz(p.I), Bvecx(p.I), Bvecy(p.I), Bvecz(p.I),
                Avec_x(p.I), Avec_y(p.I), Avec_z(p.I));
          }

          if (mask_local != 1.0) {
            // Failure inside mask
            c2p_Noble.bh_interior<EOSType, false>(eos_3p, pv_seeds, cv, glo);
            pv = pv_seeds;
          } else {
            // Failure outside, set to atmo
            cv.dBvec(0) = sqrt_detg * Bup(0);
            cv.dBvec(1) = sqrt_detg * Bup(1);
            cv.dBvec(2) = sqrt_detg * Bup(2);
            pv.Bvec = Bup;
            atmo.set(pv, cv, glo);
          }
        }
      }

      // Inside mask, C2P success
      if ((mask_local != 1.0) && c2p_flag_local) {
        c2p_Noble.bh_interior<EOSType, true>(eos_3p, pv, cv, glo);
      }
    }

    con2prim_flag(p.I) = c2p_flag_code;
    count_repair(counts, cell_atmo, set_atmo || rep_first.set_atmo ||
        rep_second.set_atmo || rep_ent.set_atmo ||
        (!c2p_flag_local && mask_local == 1.0));
    count_repair(counts, backup_failure,
        rep_second.status != c2p_report::ERR_CODE_NOT_SET && rep_second.failed());
    count_repair(counts, conservative_recompute,
        rep_first.adjust_cons || rep_second.adjust_cons || rep_ent.adjust_cons ||
        set_atmo || !c2p_flag_local);
    count_repair(counts, rho_clamp, rep_first.rho_clamped +
        rep_second.rho_clamped + rep_ent.rho_clamped);
    count_repair(counts, eps_clamp, rep_first.eps_clamped +
        rep_second.eps_clamped + rep_ent.eps_clamped);
    count_repair(counts, temp_clamp, rep_first.temp_clamped +
        rep_second.temp_clamped + rep_ent.temp_clamped);
    count_repair(counts, ye_clamp, rep_first.ye_clamped +
        rep_second.ye_clamped + rep_ent.ye_clamped);

    // ----- ----- C2P ----- -----

    // ----- Write to gfs -----

    // dummy vars
    CCTK_REAL Ex, Ey, Ez;

    // Write back pv
    pv.scatter(rho(p.I), eps(p.I), Ye(p.I), press(p.I), temperature(p.I),
               entropy(p.I), velx(p.I), vely(p.I), velz(p.I), wlor, Bvecx(p.I),
               Bvecy(p.I), Bvecz(p.I), Ex, Ey, Ez);

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
               DEnt(p.I), dBx(p.I), dBy(p.I), dBz(p.I));

    // Update saved prims
    saved_rho(p.I) = rho(p.I);
    saved_velx(p.I) = velx(p.I);
    saved_vely(p.I) = vely(p.I);
    saved_velz(p.I) = velz(p.I);
    saved_eps(p.I) = eps(p.I);
    saved_Ye(p.I) = Ye(p.I);

    // Compute diagnostics
    v_low = calc_contraction(glo, pv.vel);
    const vec<CCTK_REAL, 3> B_low = calc_contraction(glo, pv.Bvec);

    const CCTK_REAL alp_b0 = wlor * calc_contraction(pv.Bvec, v_low);

    const CCTK_REAL B2 = calc_contraction(pv.Bvec, B_low);
    const CCTK_REAL bsq = (B2 + alp_b0 * alp_b0) / (wlor * wlor);

    w_lorentz(p.I) = wlor;
    B_norm(p.I) = sqrt(B2);
    b2small(p.I) = bsq;
    volform(p.I) = sqrt_detg;
  };

  cctk_grid.loop_all_device<1, 1, 1>(grid.nghostzones, c2p_impl);
  diagnostics.report(cctkGH, "C2P");
}

extern "C" void AsterX_Con2Prim(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AsterX_Con2Prim;
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(evolution_eos, "Hybrid") && thermal_eos_atmo) {
    CCTK_ERROR("Hybrid EOS does not implement *_from_rho_temp_ye; set "
               "Con2PrimFactory::thermal_eos_atmo = no.");
  }
  if (CCTK_EQUALS(evolution_eos, "Tabulated3d") && !use_temperature) {
    CCTK_ERROR("Tabulated3d requires Con2PrimFactory::use_temperature = yes.");
  }
  if (CCTK_EQUALS(evolution_eos, "Tabulated3d") && use_press_atmo) {
    CCTK_ERROR("Tabulated3d does not support eps_from_rho_press_ye; set "
               "Con2PrimFactory::use_press_atmo = no.");
  }

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
    if (global_eos_3p_hyb_pwpoly) {
      // pwpoly cold + pwpoly-hybrid
      auto eos_cold = global_eos_1p_pwpoly;
      auto eos_3p_hyb = global_eos_3p_hyb_pwpoly;

      if (!eos_cold) {
        CCTK_ERROR("Hybrid(PWPolytropic) selected but no pwpoly cold EOS was "
                   "initialized");
      }
      AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_cold, eos_3p_hyb);

    } else if (global_eos_3p_hyb_poly) {
      // poly cold + poly-hybrid
      auto eos_cold = global_eos_1p_poly;
      auto eos_3p_hyb = global_eos_3p_hyb_poly;

      if (!eos_cold) {
        CCTK_ERROR("Hybrid(Polytropic) selected but no polytropic cold EOS was "
                   "initialized");
      }
      AsterX_Con2Prim_typeEoS(CCTK_PASS_CTOC, eos_cold, eos_3p_hyb);

    } else {
      CCTK_ERROR(
          "Hybrid EOS selected but no hybrid EOS object was initialized");
    }

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

extern "C" void AsterX_SaveC2P(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_SaveC2P;
  grid.loop_all_device<1, 1, 1>(grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) {
        rho_c2p(p.I) = rho(p.I);
        eps_c2p(p.I) = eps(p.I);
        ye_c2p(p.I) = Ye(p.I);
        vx_c2p(p.I) = velx(p.I);
        vy_c2p(p.I) = vely(p.I);
        vz_c2p(p.I) = velz(p.I);
        flag_c2p(p.I) = con2prim_flag(p.I);
      });
}

template <typename EOSIDType, typename EOSType>
void InterpolateFailed(CCTK_ARGUMENTS, const EOSIDType *eos_1p,
                       const EOSType *eos_3p) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Con2Prim_Interpolate_Failed;
  DECLARE_CCTK_PARAMETERS;
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  grid.loop_int_device<1, 1, 1>(grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) {
        if (flag_c2p(p.I) != C2P_FAIL)
          return;
        // All neighbours come from the frozen snapshot, including flags.
        const auto flag_nbs = get_neighbors(flag_c2p, p);
        const auto rho_nbs = get_neighbors(rho_c2p, p);
        const auto eps_nbs = get_neighbors(eps_c2p, p);
        const auto Ye_nbs = get_neighbors(ye_c2p, p);
        const auto velx_nbs = get_neighbors(vx_c2p, p);
        const auto vely_nbs = get_neighbors(vy_c2p, p);
        const auto velz_nbs = get_neighbors(vz_c2p, p);
        const auto good_nb = [&](int i) {
          return flag_nbs(i) != C2P_FAIL && flag_nbs(i) != C2P_INIT &&
                 std::isfinite(rho_nbs(i)) && std::isfinite(eps_nbs(i)) &&
                 std::isfinite(Ye_nbs(i)) && std::isfinite(velx_nbs(i)) &&
                 std::isfinite(vely_nbs(i)) && std::isfinite(velz_nbs(i));
        };
        CCTK_REAL sum_nbs = 0.0;
        for (int i = 0; i < 6; ++i)
          sum_nbs += good_nb(i);
        if (sum_nbs == 0.0)
          return;
        const auto average = [&](const auto &values) {
          CCTK_REAL value = 0.0;
          for (int i = 0; i < 6; ++i)
            if (good_nb(i))
              value += values(i);
          return value / sum_nbs;
        };
        const auto atmo = make_atmosphere(eos_1p, eos_3p,
            sqrt(p.x * p.x + p.y * p.y + p.z * p.z), rho_abs_min, p_atmo, t_atmo,
            Ye_atmo, r_atmo, n_rho_atmo, n_press_atmo, n_temp_atmo, atmo_tol,
            thermal_eos_atmo, use_press_atmo);
        const smat<CCTK_REAL, 3> g([&](int i, int j) ARITH_INLINE {
          return calc_avg_v2c(gf_g(i, j), p);
        });
        prim_vars pv;
        pv.Bvec = {Bvecx(p.I), Bvecy(p.I), Bvecz(p.I)};
        const CCTK_REAL rhoL = average(rho_nbs);
        if (rhoL <= atmo.rho_cut) {
          atmo.set(pv);
        } else {
          auto state = state_from_rho_eps_ye(
              eos_3p, rhoL, average(eps_nbs), average(Ye_nbs));
          if (use_press_atmo && state.press < atmo.press_atmo)
            state = state_from_rho_press_ye(
                eos_3p, state.rho, atmo.press_atmo, state.Ye);
          else if (!use_press_atmo && state.temperature < atmo.temp_atmo)
            state = state_from_rho_temp_ye(
                eos_3p, state.rho, atmo.temp_atmo, state.Ye);
          set_thermo_state(pv, state);
          pv.vel = {average(velx_nbs), average(vely_nbs), average(velz_nbs)};
          const auto v_low = calc_contraction(g, pv.vel);
          const CCTK_REAL vsq = calc_contraction(v_low, pv.vel);
          const CCTK_REAL vlim = vw_lim / sqrt(1.0 + vw_lim * vw_lim);
          if (vsq > vlim * vlim)
            pv.vel *= vlim / sqrt(vsq);
          pv.w_lor = calc_wlorentz(calc_contraction(g, pv.vel), pv.vel);
          pv.E = calc_contraction(calc_inv(g, calc_det(g)),
                                   calc_cross_product(pv.Bvec, pv.vel));
        }
        cons_vars cv;
        cv.from_prim(pv, g);
        CCTK_REAL Ex, Ey, Ez;
        pv.scatter(rho(p.I), eps(p.I), Ye(p.I), press(p.I), temperature(p.I),
            entropy(p.I), velx(p.I), vely(p.I), velz(p.I), w_lorentz(p.I),
            Bvecx(p.I), Bvecy(p.I), Bvecz(p.I), Ex, Ey, Ez);
        // The magnetic field is unchanged, so do not write the staggered B.
        dens(p.I) = cv.dens;
        momx(p.I) = cv.mom(0);
        momy(p.I) = cv.mom(1);
        momz(p.I) = cv.mom(2);
        tau(p.I) = cv.tau;
        DYe(p.I) = cv.DYe;
        DEnt(p.I) = cv.DEnt;
        saved_rho(p.I) = pv.rho;
        saved_eps(p.I) = pv.eps;
        saved_Ye(p.I) = pv.Ye;
        saved_velx(p.I) = pv.vel(0);
        saved_vely(p.I) = pv.vel(1);
        saved_velz(p.I) = pv.vel(2);
        const CCTK_REAL rhoh = pv.rho * (1.0 + pv.eps) + pv.press;
        zvec_x(p.I) = pv.w_lor * pv.vel(0);
        zvec_y(p.I) = pv.w_lor * pv.vel(1);
        zvec_z(p.I) = pv.w_lor * pv.vel(2);
        svec_x(p.I) = rhoh * pv.w_lor * pv.w_lor * pv.vel(0);
        svec_y(p.I) = rhoh * pv.w_lor * pv.w_lor * pv.vel(1);
        svec_z(p.I) = rhoh * pv.w_lor * pv.w_lor * pv.vel(2);
        const CCTK_REAL B2 = calc_contraction(pv.Bvec, calc_contraction(g, pv.Bvec));
        const CCTK_REAL Bv = calc_contraction(pv.Bvec, calc_contraction(g, pv.vel));
        B_norm(p.I) = sqrt(B2);
        b2small(p.I) = B2 / (pv.w_lor * pv.w_lor) + Bv * Bv;
        volform(p.I) = sqrt(calc_det(g));
        con2prim_flag(p.I) = C2P_AVG;
      });
}

extern "C" void AsterX_Con2Prim_Interpolate_Failed(CCTK_ARGUMENTS) {
  const auto run = [&](const auto *eos) {
    if (global_eos_1p_pwpoly)
      InterpolateFailed(cctkGH, global_eos_1p_pwpoly, eos);
    else
      InterpolateFailed(cctkGH, global_eos_1p_poly, eos);
  };
  if (global_eos_3p_ig)
    run(global_eos_3p_ig);
  else if (global_eos_3p_tab3d)
    run(global_eos_3p_tab3d);
  else
    CCTK_ERROR("Neighbour C2P repair supports IdealGas and Tabulated3d");
}

} // namespace AsterX
