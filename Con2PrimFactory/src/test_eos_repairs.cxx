#include <cctk.h>
#include <cctk_Arguments.h>

#include <array>
#include <cmath>
#include <limits>

#include "c2p_1DPalenzuela.hxx"
#include "c2p_1DRePrimAnd.hxx"
#include "c2p_2DNoble.hxx"
#include "setup_eos.hxx"

namespace Con2PrimFactory {
using namespace EOSX;

namespace {
void check_close(const char *quantity, const CCTK_REAL actual,
                  const CCTK_REAL expected, const CCTK_REAL rtol = 1.0e-7,
                  const CCTK_REAL atol = 1.0e-13) {
  if (!std::isfinite(actual) || !std::isfinite(expected) ||
      fabs(actual - expected) > atol + rtol * fabs(expected))
    CCTK_VERROR("EOS repair test: %s: actual=%.16e expected=%.16e",
                quantity, actual, expected);
}

template <typename EOSType>
void test_roundtrip(const EOSType *eos, const CCTK_REAL rho,
                     const CCTK_REAL temp, const CCTK_REAL Ye) {
  const auto state = state_from_rho_temp_ye(eos, rho, temp, Ye);
  const auto recovered = state_from_rho_eps_ye(eos, rho, state.eps, Ye);
  check_close("T -> eps -> T", recovered.temperature, state.temperature);
  check_close("eps -> T -> eps", recovered.eps, state.eps);
  check_close("round-trip P", recovered.press, state.press);
  check_close("round-trip kappa", recovered.kappa, state.kappa);
  const CCTK_REAL h = 1.0 + state.eps + state.press / state.rho;
  const auto from_h = state_from_rho_enthalpy_ye(eos, rho, h, Ye);
  if (!from_h.enthalpy_converged)
    CCTK_ERROR("EOS repair test: enthalpy inversion did not converge");
  check_close("enthalpy -> eps", from_h.state.eps, state.eps);
}

template <typename EOSType>
void test_recovery(const EOSType *eos, const thermo_state &state,
                   const atmosphere &atmo, const bool magnetized,
                   const CCTK_REAL vx = 0.1) {
  const smat<CCTK_REAL, 3> g{1.0, 0.0, 0.0, 1.0, 0.0, 1.0};
  const vec<CCTK_REAL, 3> beta{0.0, 0.0, 0.0};
  const vec<CCTK_REAL, 3> vel{vx, -0.5 * vx, 0.2 * vx};
  const CCTK_REAL B = magnetized ? sqrt(0.1 * state.press) : 0.0;
  const vec<CCTK_REAL, 3> Bvec{B, -0.5 * B, 0.25 * B};
  const CCTK_REAL wlor = 1.0 / sqrt(1.0 - calc_contraction(vel, vel));
  prim_vars seed{state.rho, state.eps, state.Ye, state.press,
                 state.temperature, state.kappa, vel, wlor, Bvec};
  seed.E = calc_cross_product(Bvec, vel);
  cons_vars original;
  original.from_prim(seed, g);

  c2p_2DNoble noble(eos, atmo, 150, 1.0e-10, 0.0, 10.0, 100.0,
                     1.0, 1.0, 10.0, 1.0e12, 1.0e12, true, false,
                     true, false, false, 1.0);
  c2p_1DPalenzuela pal(eos, atmo, 150, 1.0e-10, 0.0, 10.0, 100.0,
                       1.0, 1.0, 10.0, 1.0e12, 1.0e12, true, false,
                       true, false, false, 1.0);
  c2p_1DRePrimAnd rpa(eos, atmo, 150, 1.0e-10, 0.0, 10.0, 100.0,
                      1.0, 1.0, 10.0, 1.0e12, 1.0e12, true, false,
                      true, false, false, 1.0);
  for (int solver = 0; solver < 3; ++solver) {
    cons_vars cv = original;
    prim_vars pv;
    prim_vars trial = seed;
    trial.eps *= 1.01;
    c2p_report rep;
    if (solver == 0)
      noble.solve(eos, pv, trial, cv, 1.0, beta, g, rep);
    else if (solver == 1)
      pal.solve(eos, pv, cv, 1.0, beta, g, rep);
    else {
      auto counted = *eos;
      eos_call_counts counts;
      counted.call_counts = &counts;
      rpa.solve(&counted, pv, cv, 1.0, beta, g, rep);
      // These states need no limiter: kappa must be filled once, by the
      // finalizer, not separately before it. Check both EOS implementations.
      if (!rep.failed() &&
          counts.value[static_cast<int>(eos_call::kappa)] != 1)
        CCTK_ERROR("EOS repair test: RePrimAnd repeated thermal closure");
    }
    if (rep.failed()) {
      rep.debug_message();
      CCTK_VERROR("EOS repair test: energy C2P %d failed (magnetized=%d)",
                  solver, int(magnetized));
    }
    check_close("C2P rho", pv.rho, seed.rho);
    check_close("C2P eps", pv.eps, seed.eps);
    check_close("C2P T", pv.temperature, seed.temperature);
    check_close("C2P P", pv.press, seed.press);
    check_close("C2P kappa", pv.entropy, seed.entropy);
    check_close("C2P DEnt", cv.DEnt, cv.dens * pv.entropy);
    for (int d = 0; d < 3; ++d)
      check_close("C2P velocity", pv.vel(d), seed.vel(d));
  }

  // Compare the general-EOS Noble Jacobian with independent differences.
  const CCTK_REAL Z = (state.rho * (1.0 + state.eps) + state.press) * wlor * wlor;
  const CCTK_REAL vsq = calc_contraction(vel, vel);
  CCTK_REAL P, dPdZ, dPdvsq, Plo, Phi, unused1, unused2;
  if (!noble.get_Press_funcZVsq(P, dPdZ, dPdvsq, Z, vsq, eos, original))
    CCTK_ERROR("EOS repair test: Noble Jacobian evaluation failed");
  const CCTK_REAL dZ = 1.0e-5 * Z;
  const CCTK_REAL dv = 1.0e-6;
  if (!noble.get_Press_funcZVsq(Plo, unused1, unused2, Z-dZ, vsq, eos, original) ||
      !noble.get_Press_funcZVsq(Phi, unused1, unused2, Z+dZ, vsq, eos, original))
    CCTK_ERROR("EOS repair test: Noble Z perturbation left the EOS domain");
  check_close("Noble dP/dZ", dPdZ, (Phi-Plo)/(2*dZ), 1.0e-4);
  if (!noble.get_Press_funcZVsq(Plo, unused1, unused2, Z, vsq-dv, eos, original) ||
      !noble.get_Press_funcZVsq(Phi, unused1, unused2, Z, vsq+dv, eos, original))
    CCTK_ERROR("EOS repair test: Noble v2 perturbation left the EOS domain");
  check_close("Noble dP/dv2", dPdvsq, (Phi-Plo)/(2*dv), 1.0e-4);

  // An atmosphere reset, including magnetic energy, must be a fixed point.
  prim_vars pv = seed;
  cons_vars cv = original;
  atmo.set(pv, cv, g);
  prim_vars recovered;
  c2p_report rep;
  noble.solve(eos, recovered, pv, cv, 1.0, beta, g, rep);
  if (rep.failed() || !rep.set_atmo)
    CCTK_ERROR("EOS repair test: atmosphere was not recovered exactly");
  check_close("atmosphere eps", recovered.eps, atmo.eps_atmo, 0.0, 0.0);
  check_close("atmosphere T", recovered.temperature, atmo.temp_atmo, 0.0, 0.0);
  check_close("atmosphere DEnt", cv.DEnt, cv.dens * atmo.entropy_atmo, 0.0, 0.0);
  check_close("atmosphere tau", cv.tau,
              cv.dens * atmo.eps_atmo + 0.5 * calc_contraction(Bvec, Bvec));
}

template <typename EOSType>
void test_tau(const EOSType *eos, const thermo_state &state,
              const atmosphere &atmo) {
  // Non-unit determinant catches confusion between densitized tau and eps.
  const smat<CCTK_REAL, 3> g{2.0, 0.0, 0.0, 3.0, 0.0, 4.0};
  const CCTK_REAL sqrtg = sqrt(calc_det(g));
  prim_vars pv{state.rho, state.eps, state.Ye, state.press,
               state.temperature, state.kappa, {0.0, 0.0, 0.0}, 1.0,
               {0.001, -0.002, 0.003}};
  pv.E = {0.0, 0.0, 0.0};
  cons_vars cv;
  cv.from_prim(pv, g);
  c2p_2DNoble c2p(eos, atmo, 150, 1.0e-10, 0.0, 10.0, 100.0,
                  1.0, 1.0, 10.0, 1.0e12, 1.0e12, true, false,
                  true, false, false, 1.0);
  const CCTK_REAL valid_tau = cv.tau;
  c2p.cons_floors_and_ceilings(eos, cv, g, 1.0e-15);
  check_close("valid tau unchanged", cv.tau, valid_tau, 0.0, 0.0);
  const CCTK_REAL tau_mag = 0.5 * sqrtg *
      calc_contraction(pv.Bvec, calc_contraction(g, pv.Bvec));
  cv.tau = tau_mag + cv.dens * (fmin(0.0, eos->rgeps.min) - 0.1);
  c2p.cons_floors_and_ceilings(eos, cv, g, 1.0e-15);
  const auto er = eos->range_eps_from_rho_ye(pv.rho, pv.Ye);
  check_close("invalid tau repaired", cv.tau,
              tau_mag + cv.dens * er.min + sqrtg * 1.0e-15);
  atmo.set(pv, cv, g);
  const CCTK_REAL atmo_tau = cv.tau;
  c2p.cons_floors_and_ceilings(eos, cv, g, 1.0e-15);
  check_close("atmosphere tau unchanged", cv.tau, atmo_tau, 0.0, 0.0);
}

template <typename EOSType>
void test_temp_floor(const EOSType *eos, const thermo_state &state) {
  const auto hotter = state_from_rho_temp_ye(
      eos, state.rho, 2.0 * state.temperature, state.Ye);
  const auto heated = state_with_temp_press_floor(
      eos, state.rho, state.temperature, state.Ye, hotter.press);
  if (heated.temperature < state.temperature || heated.press < hotter.press)
    CCTK_ERROR("EOS repair test: thermal pressure floor was not reached");
  check_close("pressure floor temperature", heated.temperature, hotter.temperature);
  const auto unchanged = state_with_temp_press_floor(
      eos, state.rho, state.temperature, state.Ye, 0.5 * state.press);
  check_close("inactive pressure floor", unchanged.eps, state.eps, 0.0, 0.0);
  const auto upper = state_from_rho_temp_ye(
      eos, state.rho, eos->rgtemp.max, state.Ye);
  const auto saturated = state_with_temp_press_floor(
      eos, state.rho, state.temperature, state.Ye, 2.0 * upper.press);
  check_close("unreachable floor saturates at Tmax", saturated.temperature,
              eos->rgtemp.max, 0.0, 0.0);
}

// Exercise the shared finalizer without depending on root convergence.
struct test_c2p : c2p_2DNoble {
  using c2p_2DNoble::c2p_2DNoble;
  using c2p::prims_floors_and_ceilings;
};

template <typename EOSType>
void test_mag_floors(const EOSType *eos, const atmosphere &atmo,
                     const bool use_temp) {
  const smat<CCTK_REAL, 3> g{2.0, 0.0, 0.0, 3.0, 0.0, 4.0};
  const vec<CCTK_REAL, 3> beta{0.0, 0.0, 0.0};
  const auto state = state_from_rho_temp_ye(
      eos, 0.25 * eos->rgrho.max, 0.25 * eos->rgtemp.max,
      0.5 * (eos->rgye.min + eos->rgye.max));
  const auto upper = state_from_rho_temp_ye(
      eos, state.rho, eos->rgtemp.max, state.Ye);
  const CCTK_REAL eps = std::numeric_limits<CCTK_REAL>::epsilon();
  const CCTK_REAL tol = 64.0 * eps;

  auto check = [&](const char *name, const CCTK_REAL rho_min,
                   const CCTK_REAL press_min, const vec<CCTK_REAL, 3> &vel,
                   const bool failed, const bool magnetized = true) {
    const vec<CCTK_REAL, 3> Bvec{
        magnetized ? 1.0 / sqrt(g(0, 0)) : 0.0, 0.0, 0.0};
    const auto v_low = calc_contraction(g, vel);
    const CCTK_REAL vsq = calc_contraction(vel, v_low);
    const CCTK_REAL B2 =
        calc_contraction(Bvec, calc_contraction(g, Bvec));
    const CCTK_REAL Bdotv = calc_contraction(Bvec, v_low);
    const CCTK_REAL bsq = B2 * (1.0 - vsq) + Bdotv * Bdotv;
    // Choose the limits to request the supplied density and pressure floors.
    const CCTK_REAL sigma = (magnetized ? bsq : 1.0) / rho_min;
    const CCTK_REAL inv_beta = 0.5 * (magnetized ? bsq : 1.0) / press_min;
    prim_vars pv{state.rho, state.eps, state.Ye, state.press,
                 state.temperature, state.kappa, vel, 1.0 / sqrt(1.0 - vsq),
                 Bvec};
    pv.E = calc_contraction(calc_inv(g, calc_det(g)),
                            calc_cross_product(Bvec, vel));
    cons_vars cv;
    cv.from_prim(pv, g);
    test_c2p c2p(eos, atmo, 150, 1.0e-10, 0.0, 10.0, 100.0,
                  1.0, 1.0, 10.0, sigma, inv_beta, true, false,
                  use_temp, false, false, 1.0);
    c2p_report rep;
    rep.status = c2p_report::SUCCESS;
    c2p.prims_floors_and_ceilings(eos, pv, cv, 1.0, beta, g, rep);
    if (failed) {
      if (rep.status != c2p_report::B_LIMIT)
        CCTK_VERROR("EOS repair test: %s did not report B_LIMIT (use_temp=%d)",
                    name, int(use_temp));
      return;
    }
    if (rep.failed()) {
      rep.debug_message();
      CCTK_VERROR("EOS repair test: %s failed (use_temp=%d)",
                  name, int(use_temp));
    }

    const auto er = eos->range_eps_from_rho_ye(pv.rho, pv.Ye);
    const CCTK_REAL Bdotv_new =
        calc_contraction(Bvec, calc_contraction(g, pv.vel));
    const CCTK_REAL vsq_new =
        calc_contraction(pv.vel, calc_contraction(g, pv.vel));
    const CCTK_REAL bsq_new =
        B2 * (1.0 - vsq_new) + Bdotv_new * Bdotv_new;
    if (!std::isfinite(pv.rho) || !std::isfinite(pv.eps) ||
        !std::isfinite(pv.press) || !std::isfinite(bsq_new) ||
        pv.rho < eos->rgrho.min || pv.rho > eos->rgrho.max ||
        pv.eps < er.min || pv.eps > er.max ||
        pv.rho < (1.0 - tol) * rho_min ||
        pv.press < (1.0 - tol) * press_min ||
        bsq_new > (1.0 + tol) * sigma * pv.rho ||
        bsq_new > (1.0 + tol) * 2.0 * inv_beta * pv.press)
      CCTK_VERROR("EOS repair test: %s lost an EOS or magnetic limit", name);
    check_close("magnetic floor Lorentz factor", pv.w_lor,
                1.0 / sqrt(1.0 - vsq_new));
    const auto closed = state_from_rho_eps_ye(eos, pv.rho, pv.eps, pv.Ye);
    check_close("magnetic floor pressure", pv.press, closed.press);
    check_close("magnetic floor temperature", pv.temperature, closed.temperature);
    check_close("magnetic floor kappa", pv.entropy, closed.kappa);
    for (int d = 0; d < 3; ++d)
      check_close("magnetic floor B unchanged", pv.Bvec(d), Bvec(d), 0.0, 0.0);
  };

  const std::array<vec<CCTK_REAL, 3>, 4> velocities{
      vec<CCTK_REAL, 3>{0.0, 0.0, 0.0},
      vec<CCTK_REAL, 3>{0.3, 0.0, 0.0},
      vec<CCTK_REAL, 3>{0.0, 0.2, 0.0},
      vec<CCTK_REAL, 3>{0.3, 0.2, 0.1}};
  for (const auto &vel : velocities) {
    check("zero field", state.rho, state.press, vel, false, false);
    check("inactive floors", 0.5 * state.rho, 0.5 * state.press, vel, false);
    check("density floor", 2.0 * state.rho, state.press, vel, false);
    check("pressure floor", state.rho, 2.0 * state.press, vel, false);
    check("both floors", 2.0 * state.rho, 2.0 * state.press, vel, false);
    check("density saturation", 2.0 * eos->rgrho.max, state.press, vel, true);
    check("pressure saturation", state.rho, 2.0 * upper.press, vel, true);
    check("density roundoff", (1.0 + 8.0 * eps) * eos->rgrho.max,
          state.press, vel, false);
    check("pressure roundoff", state.rho,
          (1.0 + 8.0 * eps) * upper.press, vel, false);
    check("density beyond roundoff", (1.0 + 1.0e-8) * eos->rgrho.max,
          state.press, vel, true);
    check("pressure beyond roundoff", state.rho,
          (1.0 + 1.0e-8) * upper.press, vel, true);
  }
}

void test_ideal_gas() {
  eos_3p_idealgas eos;
  eos_3p::range rgeps{0.0, 2.0}, rgrho{1.0e-14, 1.0}, rgye{0.0, 1.0};
  eos.init(2.0, 1.0, rgeps, rgrho, rgye);
  eos_1p_polytropic cold;
  cold.init(2.0, 100.0, 1.0);
  const auto inner = make_atmosphere(&cold, &eos, 1.0, 1.0e-6, 1.0e-8,
      1.0e-4, 0.5, 1.0, 6.0, 12.0, 6.0, 1.0e-3, false, false);
  const auto outer = make_atmosphere(&cold, &eos, 2.0, 1.0e-6, 1.0e-8,
      1.0e-4, 0.5, 1.0, 6.0, 12.0, 6.0, 1.0e-3, false, false);
  check_close("cold density grading", outer.rho_atmo, inner.rho_atmo / 64.0);
  check_close("cold energy grading", outer.eps_atmo, inner.eps_atmo / 64.0);
  // The outer pressure is below the default absolute tolerance.
  check_close("cold pressure grading", outer.press_atmo,
              inner.press_atmo / 4096.0, 1.0e-12, 0.0);
  check_close("cold temperature grading", outer.temp_atmo, inner.temp_atmo / 64.0);
  const auto state = state_from_rho_temp_ye(&eos, 1.0e-3, 0.02, 0.5);
  const auto pressure_state = state_from_rho_press_ye(&eos, state.rho,
                                                     state.press, state.Ye);
  check_close("ideal-gas pressure closure", pressure_state.eps, state.eps);
  check_close("ideal-gas kappa", state.kappa, state.press / (state.rho*state.rho));
  test_roundtrip(&eos, state.rho, state.temperature, state.Ye);
  test_recovery(&eos, state, inner, false);
  test_recovery(&eos, state, inner, true);
  test_tau(&eos, state, inner);
  test_temp_floor(&eos, state);
  for (const CCTK_REAL gamma : {1.4, 5.0 / 3.0, 2.0}) {
    eos.init(gamma, 1.0, rgeps, rgrho, rgye);
    const auto atmo = make_atmosphere(&cold, &eos, 1.0, 1.0e-6, 0.0,
        1.0e-4, 0.5, 1.0, 0.0, 0.0, 0.0, 1.0e-3, true, false);
    test_mag_floors(&eos, atmo, true);
    test_mag_floors(&eos, atmo, false);
    for (const CCTK_REAL rho : {1.0e-4, 1.0e-2})
      for (const CCTK_REAL vx : {0.1, 0.5, 0.7}) {
        const auto sample = state_from_rho_temp_ye(&eos, rho, 0.02, 0.5);
        test_recovery(&eos, sample, atmo, false, vx);
        test_recovery(&eos, sample, atmo, true, vx);
      }
  }
}

void test_shifted_table() {
  // Exact log-linear synthetic EOS: eps=T-shift, P=rho*T. No derivative
  // columns are supplied, as in the current CompOSE reader.
  std::array<CCTK_REAL, 3> lr{log(1.0e-6), log(1.0e-4), log(1.0e-2)};
  std::array<CCTK_REAL, 3> lt{log(1.0e-3), log(1.0e-2), log(1.0e-1)};
  std::array<CCTK_REAL, 2> ye{0.1, 0.5};
  std::array<CCTK_REAL, 18 * NTABLES> data{};
  CCTK_REAL shift = 0.05;
  for (int k = 0; k < 2; ++k)
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i) {
        const int offset = NTABLES * (i + 3*(j + 3*k));
        data[offset + eos_3p_tabulated3d::PRESS] = lr[i] + lt[j];
        data[offset + eos_3p_tabulated3d::EPS] = lt[j];
        data[offset + eos_3p_tabulated3d::S] = lt[j] - lr[i];
        const CCTK_REAL temp = exp(lt[j]);
        data[offset + eos_3p_tabulated3d::CS2] = 2*temp/(1-shift+2*temp);
      }
  linear_interp_uniform_ND_t<CCTK_REAL, 3, NTABLES> interp(
      data.data(), {3, 3, 2}, lr.data(), lt.data(), ye.data());
  eos_3p_tabulated3d eos;
  eos.interptable = &interp;
  eos.energy_shift = &shift;
  eos.rgrho = {exp(lr.front()), exp(lr.back())};
  eos.rgtemp = {exp(lt.front()), exp(lt.back())};
  eos.rgye = {ye.front(), ye.back()};
  eos.rgeps = eos.compute_eps_range_full_table();
  check_close("physical table minimum", eos.rgeps.min, 0.001-shift);
  const auto state = state_from_rho_temp_ye(&eos, 1.0e-3, 0.02, 0.3);
  check_close("negative physical eps", state.eps, -0.03);
  CCTK_REAL P, dpdrho, dpdeps;
  eos.press_derivs_from_rho_eps_ye(P, dpdrho, dpdeps, state.rho, state.eps, state.Ye);
  check_close("table dP/drho at fixed eps", dpdrho, state.temperature);
  check_close("table dP/deps at fixed rho", dpdeps, state.rho);
  test_roundtrip(&eos, state.rho, state.temperature, state.Ye);
  test_roundtrip(&eos, state.rho, eos.rgtemp.min, state.Ye);
  test_roundtrip(&eos, state.rho, eos.rgtemp.max, state.Ye);
  eos_1p_polytropic cold;
  cold.init(2.0, 100.0, 1.0);
  const auto inner = make_atmosphere(&cold, &eos, 1.0, 1.0e-4, 0.0,
      0.004, 0.3, 1.0, 1.0, 0.0, 1.0, 1.0e-3, true, false);
  const auto outer = make_atmosphere(&cold, &eos, 2.0, 1.0e-4, 0.0,
      0.004, 0.3, 1.0, 1.0, 0.0, 1.0, 1.0e-3, true, false);
  check_close("tabulated density grading", outer.rho_atmo, inner.rho_atmo/2);
  check_close("tabulated temperature grading", outer.temp_atmo, inner.temp_atmo/2);
  check_close("tabulated graded energy", outer.eps_atmo, outer.temp_atmo-shift);
  test_recovery(&eos, state, inner, false);
  test_recovery(&eos, state, inner, true);
  test_tau(&eos, state, inner);
  test_temp_floor(&eos, state);
  test_mag_floors(&eos, inner, true);
  for (const CCTK_REAL rho : {3.0e-4, 3.0e-3})
    for (const CCTK_REAL temp : {0.006, 0.02, 0.07})
      for (const CCTK_REAL vx : {0.1, 0.5, 0.7}) {
        const auto sample = state_from_rho_temp_ye(&eos, rho, temp, 0.3);
        test_recovery(&eos, sample, inner, false, vx);
        test_recovery(&eos, sample, inner, true, vx);
      }
}
} // namespace

extern "C" void Con2PrimFactory_TestEOSRepairs(CCTK_ARGUMENTS) {
  test_ideal_gas();
  test_shifted_table();
  if (global_eos_3p_tab3d) {
    const auto *eos = global_eos_3p_tab3d;
    for (const CCTK_REAL f : {0.25, 0.5, 0.75}) {
      const CCTK_REAL rho = exp((1-f)*log(eos->rgrho.min) + f*log(eos->rgrho.max));
      const CCTK_REAL temp = exp((1-f)*log(eos->rgtemp.min) + f*log(eos->rgtemp.max));
      const CCTK_REAL Ye = (1-f)*eos->rgye.min + f*eos->rgye.max;
      test_roundtrip(eos, rho, temp, Ye);
    }
  }
  CCTK_INFO("EOS repair closure, grading, Jacobian and recovery tests passed");
}
} // namespace Con2PrimFactory
