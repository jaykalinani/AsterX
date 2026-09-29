#include <cctk.h>
#include <cctk_Arguments.h>
#include <AMReX_GpuLaunch.H>
#include <AMReX_GpuMemory.H>

#include "setup_eos.hxx"
#include "thermo_state.hxx"

namespace EOSX {
namespace {
void check(const char *name, const CCTK_REAL a, const CCTK_REAL b,
           const CCTK_REAL tol = 1.0e-10) {
  if (!std::isfinite(a) || !std::isfinite(b) ||
      fabs(a-b) > tol * fmax(1.0, fabs(b)))
    CCTK_VERROR("EOSX test %s failed: %.16e != %.16e", name, a, b);
}

template <bool uniform> void test_interp() {
  using interp_t = lintp_ND_t<CCTK_REAL, 3, 1, uniform>;
  std::array<CCTK_REAL, 3> x{0.0, uniform ? 0.5 : 0.25, 1.0};
  std::array<CCTK_REAL, 3> y{-1.0, uniform ? 0.0 : -0.3, 1.0};
  std::array<CCTK_REAL, 3> z{1.0, uniform ? 2.0 : 1.4, 3.0};
  std::array<CCTK_REAL, 27> data;
  const auto f = [](CCTK_REAL a, CCTK_REAL b, CCTK_REAL c) {
    return 2 + 3*a - 4*b + 5*c + 1.3*a*b + 0.7*a*b*c;
  };
  for (int k = 0; k < 3; ++k)
    for (int j = 0; j < 3; ++j)
      for (int i = 0; i < 3; ++i)
        data[i + 3*(j + 3*k)] = f(x[i], y[j], z[k]);
  interp_t interp(data.data(), {3, 3, 3}, x.data(), y.data(), z.data());
  for (const CCTK_REAL a : {0.0, 0.19, 0.5, 0.81, 1.0})
    for (const CCTK_REAL b : {-1.0, -0.37, 0.21, 1.0})
      for (const CCTK_REAL c : {1.0, 1.37, 2.43, 3.0}) {
        const auto v = interp.template interpolate_with_derivs<0>(a, b, c);
        check("interpolated value", v[0], f(a, b, c));
        check("d/dx", v[1], 3 + 1.3*b + 0.7*b*c);
        check("d/dy", v[2], -4 + 1.3*a + 0.7*a*c);
        check("d/dz", v[3], 5 + 0.7*a*b);
      }
}

void test_ideal() {
  for (const CCTK_REAL gamma : {1.4, 5.0/3.0, 2.0}) {
    eos_3p_idealgas eos;
    eos_3p::range er{0.0, 3.0}, rr{1.0e-12, 1.0}, yr{0.0, 1.0};
    eos.init(gamma, 2.0, er, rr, yr);
    eos_call_counts counts;
    eos.call_counts = &counts;
    for (const CCTK_REAL rho : {1.0e-10, 1.0e-4, 0.5})
      for (const CCTK_REAL eps : {0.0, 1.0e-6, 0.1, 2.9}) {
        const auto state = state_from_rho_eps_ye(&eos, rho, eps, 0.5);
        check("ideal pressure", state.press, (gamma-1)*rho*eps);
        check("ideal temperature", state.temperature, 2*(gamma-1)*eps);
        check("ideal kappa", state.kappa, (gamma-1)*eps*pow(rho, 1-gamma));
        check("ideal sound speed", state.cs2, gamma*(gamma-1)*eps/(1+gamma*eps));
        const auto from_t = state_from_rho_temp_ye(&eos, rho, state.temperature, 0.5);
        const auto from_p = state_from_rho_press_ye(&eos, rho, state.press, 0.5);
        check("ideal T round trip", from_t.eps, eps);
        check("ideal P round trip", from_p.eps, eps);
        const auto from_h = state_from_rho_enthalpy_ye(&eos, rho, 1+gamma*eps, 0.5);
        if (!from_h.enthalpy_converged)
          CCTK_ERROR("EOSX test: ideal enthalpy inversion failed");
        check("ideal h round trip", from_h.state.eps, eps);
      }
    const auto limited = state_from_rho_temp_ye(&eos, -1.0, -1.0, -1.0);
    check("rho lower bound", limited.rho, eos.rgrho.min, 0.0);
    check("T lower bound", limited.temperature, eos.rgtemp.min, 0.0);
    check("Ye lower bound", limited.Ye, eos.rgye.min, 0.0);
    const auto upper = state_from_rho_eps_ye(&eos, 2.0, 4.0, 2.0);
    check("rho upper bound", upper.rho, eos.rgrho.max, 0.0);
    check("eps upper bound", upper.eps, eos.rgeps.max, 0.0);
    check("Ye upper bound", upper.Ye, eos.rgye.max, 0.0);
    const auto failed = state_from_rho_enthalpy_ye(
        &eos, 0.01, std::numeric_limits<CCTK_REAL>::quiet_NaN(), 0.5);
    if (failed.enthalpy_converged)
      CCTK_ERROR("EOSX test: non-finite enthalpy was accepted");
    if (!counts.value[static_cast<int>(eos_call::pressure)] ||
        !counts.value[static_cast<int>(eos_call::enthalpy)] ||
        !counts.value[static_cast<int>(eos_call::derivatives)])
      CCTK_ERROR("EOSX test: EOS API counters did not record calls");
  }
}

void test_table_closure() {
  // A shifted log-linear table with rho- and Ye-dependent energy bounds.
  // This exercises the real interpolator and inverse, not a mock EOS.
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
        data[offset + eos_3p_tabulated3d::EPS] = lt[j] + 0.1*lr[i] + 0.2*ye[k];
        data[offset + eos_3p_tabulated3d::S] = lt[j] - lr[i] + ye[k];
        data[offset + eos_3p_tabulated3d::CS2] = 0.2;
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
  eos_call_counts counts;
  eos.call_counts = &counts;

  const auto check_value = [](const char *name, CCTK_REAL a, CCTK_REAL b) {
    if (!std::isfinite(a) || !std::isfinite(b) ||
        fabs(a-b) > 1.0e-10*fmax(1.0e-12, fabs(b)))
      CCTK_VERROR("EOSX test %s failed: %.16e != %.16e", name, a, b);
  };
  for (const CCTK_REAL rho : {1.0e-8, 1.0e-6, 3.0e-4, 1.0e-2, 1.0})
    for (const CCTK_REAL temp : {0.0, 0.001, 0.007, 0.1, 1.0})
      for (const CCTK_REAL Ye : {-0.1, 0.1, 0.3, 0.5, 0.9}) {
        counts = {};
        const auto state = state_from_rho_temp_ye(&eos, rho, temp, Ye);
        if (counts.value[static_cast<int>(eos_call::table_inverse)] != 0)
          CCTK_ERROR("EOSX test: temperature closure inverted the table");
        const CCTK_REAL r = limit_to_range(rho, eos.rgrho);
        const CCTK_REAL t = limit_to_range(temp, eos.rgtemp);
        const CCTK_REAL y = limit_to_range(Ye, eos.rgye);
        check_value("table T-primary eps", state.eps,
                    t*pow(r, 0.1)*exp(0.2*y)-shift);
        check_value("table T-primary P", state.press, r*t);
        check_value("table T-primary kappa", state.kappa, log(t)-log(r)+y);
        check_value("table T-primary cs2", state.cs2, 0.2);
        check("table T authority", state.temperature, t, 0.0);

        const auto er = eos.range_eps_from_rho_ye(r, y);
        for (const CCTK_REAL eps : {er.min-1.0, er.min, state.eps,
                                    er.max, er.max+1.0}) {
          counts = {};
          const auto recovered = state_from_rho_eps_ye(&eos, rho, eps, Ye);
          if (counts.value[static_cast<int>(eos_call::table_inverse)] != 1)
            CCTK_ERROR("EOSX test: energy closure did not invert exactly once");
          CCTK_REAL eps_ref = limit_to_range(eps, er);
          // Independent legacy calls provide a reference for each output.
          const CCTK_REAL press = eos.press_from_rho_eps_ye(r, eps_ref, y);
          const CCTK_REAL temperature = eos.temp_from_rho_eps_ye(r, eps_ref, y);
          const CCTK_REAL kappa = eos.kappa_from_rho_eps_ye(r, eps_ref, y);
          const CCTK_REAL cs = eos.csnd_from_rho_eps_ye(r, eps_ref, y);
          check_value("table energy authority", recovered.eps, eps_ref);
          check_value("table reused P", recovered.press, press);
          check_value("table reused T", recovered.temperature, temperature);
          check_value("table reused kappa", recovered.kappa, kappa);
          check_value("table reused cs2", recovered.cs2, cs*cs);
        }
      }
}

template <typename EOSType> void test_device(const EOSType *eos) {
  amrex::Gpu::DeviceScalar<unsigned int> failures(0);
  auto *failed = failures.dataPtr();
  // The active EOS lives in managed storage and is valid on the device.
  amrex::ParallelFor(64, [=] AMREX_GPU_DEVICE(int i) {
    const CCTK_REAL f = (i + 0.5) / 64.0;
    const CCTK_REAL rho = exp((1-f)*log(eos->rgrho.min) + f*log(eos->rgrho.max));
    const CCTK_REAL temp = eos->rgtemp.min + f*(eos->rgtemp.max-eos->rgtemp.min);
    const CCTK_REAL Ye = eos->rgye.min + f*(eos->rgye.max-eos->rgye.min);
    const auto state = state_from_rho_temp_ye(eos, rho, temp, Ye);
    const auto recovered = state_from_rho_eps_ye(eos, rho, state.eps, Ye);
    const CCTK_REAL actual[] = {recovered.eps, recovered.press,
                                recovered.kappa, recovered.cs2};
    const CCTK_REAL expected[] = {state.eps, state.press, state.kappa, state.cs2};
    for (int j = 0; j < 4; ++j)
      if (!std::isfinite(actual[j]) || !std::isfinite(expected[j]) ||
          fabs(actual[j]-expected[j]) > 1.0e-7*fmax(fabs(expected[j]), 1.0e-20))
        amrex::HostDevice::Atomic::Add(failed, 1U);
    if (!std::isfinite(recovered.temperature) ||
        fabs(recovered.temperature-temp) > 1.0e-7*fmax(temp, 1.0e-12))
      amrex::HostDevice::Atomic::Add(failed, 1U);
  });
  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("EOSX test: active-EOS device round trips failed");

  if constexpr (EOSType::temperature_primary) {
    amrex::Gpu::DeviceScalar<eos_call_counts> calls(eos_call_counts{});
    auto *counts = calls.dataPtr();
    amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE(int) {
      auto local = *eos;
      local.call_counts = counts;
      const CCTK_REAL rho = sqrt(eos->rgrho.min * eos->rgrho.max);
      const CCTK_REAL temp = sqrt(eos->rgtemp.min * eos->rgtemp.max);
      const CCTK_REAL Ye = 0.5 * (eos->rgye.min + eos->rgye.max);
      const auto state = state_from_rho_temp_ye(&local, rho, temp, Ye);
      if (counts->value[static_cast<int>(eos_call::table_inverse)] != 0)
        amrex::HostDevice::Atomic::Add(failed, 1U);
      const auto recovered = state_from_rho_eps_ye(&local, rho, state.eps, Ye);
      if (!std::isfinite(recovered.temperature))
        amrex::HostDevice::Atomic::Add(failed, 1U);
    });
    amrex::Gpu::streamSynchronize();
    if (failures.dataValue() ||
        calls.dataValue().value[static_cast<int>(eos_call::table_inverse)] != 1)
      CCTK_ERROR("EOSX test: device closure did not reuse table temperature");
  }
}
} // namespace

extern "C" void EOSX_Test(CCTK_ARGUMENTS) {
  test_interp<true>();
  test_interp<false>();
  test_ideal();
  test_table_closure();
  if (global_eos_3p_ig)
    test_device(global_eos_3p_ig);
  if (global_eos_3p_tab3d)
    test_device(global_eos_3p_tab3d);
  CCTK_INFO("EOSX interpolation, closure, counters and device tests passed");
}
} // namespace EOSX
