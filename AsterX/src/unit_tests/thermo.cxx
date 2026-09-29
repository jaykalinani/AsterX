#include <cctk.h>

#include "atmo.hxx"
#include "setup_eos.hxx"
#include "../eigenvalues.hxx"
#include "../atmo_cache.hxx"
#include "../test.hxx"

void AsterXTests::test_thermo() {
  using namespace EOSX;
  using namespace Con2PrimFactory;
  eos_3p_idealgas eos;
  eos_3p::range er{0.0, 2.0}, rr{1.0e-14, 1.0}, yr{0.0, 1.0};
  eos.init(2.0, 1.0, er, rr, yr);
  eos_1p_polytropic cold;
  cold.init(2.0, 100.0, 1.0);

  // A face at r=1.5, dx=1 has neighbours at r=1 and r=2. Each side
  // must use its own atmosphere, not the atmosphere evaluated at the face.
  const atmosphere atmo[2] = {
      make_atmosphere(&cold, &eos, 1.0, 1.0e-6, 0.0, 0.0, 0.5,
          1.0, 6.0, 0.0, 0.0, 1.0e-3, false, false),
      make_atmosphere(&cold, &eos, 2.0, 1.0e-6, 0.0, 0.0, 0.5,
          1.0, 6.0, 0.0, 0.0, 1.0e-3, false, false)};
  if (!isapprox(atmo[0].rho_atmo / atmo[1].rho_atmo, 64.0) ||
      !isapprox(atmo[0].press_atmo / atmo[1].press_atmo, 4096.0))
    CCTK_ERROR("AsterX test: left/right atmosphere grading failed");

  vec<CCTK_REAL, 2> rho, cs2, h;
  for (int f = 0; f < 2; ++f) {
    const auto state = state_from_rho_press_ye(
        &eos, atmo[f].rho_atmo, atmo[f].press_atmo, atmo[f].ye_atmo);
    if (!isapprox(state.eps, atmo[f].eps_atmo) ||
        !isapprox(state.kappa, atmo[f].entropy_atmo))
      CCTK_ERROR("AsterX test: pressure face closure failed");
    rho(f) = state.rho;
    cs2(f) = state.cs2;
    h(f) = 1.0 + state.eps + state.press / state.rho;
  }
  // At rest, without B, the two fast speeds reduce to +/-cs on each side.
  const auto lambda = AsterX::eigenvalues(
      1.0, 0.0, 1.0, {0.0, 0.0}, rho, cs2, {1.0, 1.0}, h, {0.0, 0.0});
  for (int f = 0; f < 2; ++f)
    if (!isapprox(lambda(f)(0), sqrt(cs2(f))) ||
        !isapprox(lambda(f)(2), -sqrt(cs2(f))))
      CCTK_ERROR("AsterX test: face-state sound speeds failed");

  // Exercise cold, temperature-primary and pressure-primary cache entries.
  // The comparison is exact: a cache hit must only copy the same state.
  for (int mode = 0; mode < 3; ++mode) {
    const bool thermal = mode != 0;
    const bool use_press = mode == 2;
    const auto inner = make_atmosphere(&cold, &eos, 0.0, 1.0e-6,
        1.0e-10, 1.0e-4, 0.5, 1.0, 0.0, 0.0, 0.0, 1.0e-3,
        thermal, use_press);
    AsterX::atmo_cache cache;
    atmosphere cached{};
    const auto load = [&](CCTK_REAL nr, CCTK_REAL np, CCTK_REAL nt,
                           CCTK_REAL tol, bool th, bool pr) {
      return cache.load(&cold, &eos, 1.0e-6, 1.0e-10, 1.0e-4, 0.5,
                         nr, np, nt, tol, th, pr, cached);
    };
    if (load(0.0, 0.0, 0.0, 1.0e-3, thermal, use_press))
      CCTK_ERROR("AsterX test: uninitialized atmosphere cache hit");
    cache.store(&cold, &eos, inner, 1.0e-6, 1.0e-10, 1.0e-4, 0.5,
                  thermal, use_press);
    for (const CCTK_REAL r : {0.0, 0.5, 1.0, 2.0, 100.0}) {
      const auto direct = make_atmosphere(&cold, &eos, r, 1.0e-6,
          1.0e-10, 1.0e-4, 0.5, 1.0, 0.0, 0.0, 0.0, 1.0e-3,
          thermal, use_press);
      if (!load(0.0, 0.0, 0.0, 1.0e-3, thermal, use_press) ||
          cached.rho_atmo != direct.rho_atmo ||
          cached.eps_atmo != direct.eps_atmo ||
          cached.ye_atmo != direct.ye_atmo ||
          cached.press_atmo != direct.press_atmo ||
          cached.temp_atmo != direct.temp_atmo ||
          cached.entropy_atmo != direct.entropy_atmo ||
          cached.rho_cut != direct.rho_cut)
        CCTK_ERROR("AsterX test: cached atmosphere differs from builder");
    }
    if (!load(0.0, 0.0, 0.0, 0.02, thermal, use_press) ||
        cached.rho_cut != inner.rho_atmo * 1.02)
      CCTK_ERROR("AsterX test: cached atmosphere ignored steered cutoff");
    if (load(1.0, 0.0, 0.0, 1.0e-3, thermal, use_press) ||
        load(0.0, 0.0, 0.0, 1.0e-3, !thermal, use_press) ||
        load(0.0, 0.0, 0.0, 1.0e-3, thermal, !use_press))
      CCTK_ERROR("AsterX test: cache ignored density grading or a mode change");
    if (thermal && load(0.0, use_press ? 1.0 : 0.0,
                          use_press ? 0.0 : 1.0, 1.0e-3, thermal, use_press))
      CCTK_ERROR("AsterX test: cache ignored active thermal grading");
    if (!load(0.0, use_press && thermal ? 0.0 : 1.0,
                !use_press && thermal ? 0.0 : 1.0,
                1.0e-3, thermal, use_press))
      CCTK_ERROR("AsterX test: inactive grading prevented a cache hit");
    auto other = eos;
    if (cache.load(&cold, &other, 1.0e-6, 1.0e-10, 1.0e-4, 0.5,
                     0.0, 0.0, 0.0, 1.0e-3, thermal, use_press, cached) ||
        cache.load(&cold, &eos, 2.0e-6, 1.0e-10, 1.0e-4, 0.5,
                     0.0, 0.0, 0.0, 1.0e-3, thermal, use_press, cached) ||
        cache.load(&cold, &eos, 1.0e-6, 1.0e-10, 2.0e-4, 0.5,
                     0.0, 0.0, 0.0, 1.0e-3, thermal, use_press, cached))
      CCTK_ERROR("AsterX test: cache accepted a changed EOS or atmosphere");
  }
  CCTK_INFO("AsterX face closure, grading and sound-speed tests passed");
}
