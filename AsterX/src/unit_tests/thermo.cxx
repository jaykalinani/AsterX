#include <cctk.h>

#include "atmo.hxx"
#include "setup_eos.hxx"
#include "../eigenvalues.hxx"
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
  CCTK_INFO("AsterX face closure, grading and sound-speed tests passed");
}
