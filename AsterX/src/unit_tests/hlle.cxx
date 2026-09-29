#include <cctk.h>

#include "../fluxes.hxx"
#include "../test.hxx"

void AsterXTests::test_hlle() {
  using namespace AsterX;
  const vec<CCTK_REAL, 2> var{1.0, 3.0};
  const vec<CCTK_REAL, 2> flux{2.0, 4.0};
  vec<vec<CCTK_REAL, 4>, 2> lam;
  for (int f = 0; f < 2; ++f)
    for (int n = 0; n < 4; ++n)
      lam(f)(n) = 0.0;
  if (hlle(lam, var, flux) != 3.0)
    CCTK_ERROR("HLLE zero-speed limit failed");
  for (int f = 0; f < 2; ++f)
    for (int n = 0; n < 4; ++n)
      lam(f)(n) = 1.0;
  if (hlle(lam, var, flux) != flux(0))
    CCTK_ERROR("HLLE right-going upwind limit failed");
  for (int f = 0; f < 2; ++f)
    for (int n = 0; n < 4; ++n)
      lam(f)(n) = -1.0;
  if (hlle(lam, var, flux) != flux(1))
    CCTK_ERROR("HLLE left-going upwind limit failed");
  lam(0)(0) = 1.0;
  if (hlle(lam, var, flux) != laxf(lam, var, flux))
    CCTK_ERROR("HLLE symmetric-speed limit failed");
  CCTK_INFO("HLLE limiting cases passed");
}
