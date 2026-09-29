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
  if (!std::isfinite(a) || fabs(a-b) > tol * fmax(1.0, fabs(b)))
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
    if (!std::isfinite(recovered.temperature) ||
        fabs(recovered.temperature-temp) > 1.0e-7*fmax(temp, 1.0e-12))
      amrex::HostDevice::Atomic::Add(failed, 1U);
  });
  amrex::Gpu::streamSynchronize();
  if (failures.dataValue())
    CCTK_ERROR("EOSX test: active-EOS device round trips failed");
}
} // namespace

extern "C" void EOSX_Test(CCTK_ARGUMENTS) {
  test_interp<true>();
  test_interp<false>();
  test_ideal();
  if (global_eos_3p_ig)
    test_device(global_eos_3p_ig);
  if (global_eos_3p_tab3d)
    test_device(global_eos_3p_tab3d);
  CCTK_INFO("EOSX interpolation, closure, counters and device tests passed");
}
} // namespace EOSX
