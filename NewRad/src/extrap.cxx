#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace NewRad {
using namespace Loop;

void extrap(const cGH *restrict const cctkGH, CCTK_REAL *restrict const var) {
  const GF3D<CCTK_REAL, 0, 0, 0> var_(cctkGH, var);

  loop_bnd<0, 0, 0>(cctkGH,
                    [&](const PointDesc &p) { var_(p.I) = var_(p.I - p.NI); });
}

extern "C" void NewRad_Extrapolate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_NewRad_Extrapolate;

  extrap(cctkGH, lambdaU0GF);
  extrap(cctkGH, lambdaU1GF);
  extrap(cctkGH, lambdaU2GF);

  // extrap(cctkGH, aDD00GF);
  // extrap(cctkGH, aDD01GF);
  // extrap(cctkGH, aDD02GF);
  // extrap(cctkGH, aDD11GF);
  // extrap(cctkGH, aDD12GF);
  // extrap(cctkGH, aDD22GF);
  // extrap(cctkGH, alphaGF);
  // extrap(cctkGH, betU0GF);
  // extrap(cctkGH, betU1GF);
  // extrap(cctkGH, betU2GF);
  // extrap(cctkGH, cfGF);
  // extrap(cctkGH, hDD00GF);
  // extrap(cctkGH, hDD01GF);
  // extrap(cctkGH, hDD02GF);
  // extrap(cctkGH, hDD11GF);
  // extrap(cctkGH, hDD12GF);
  // extrap(cctkGH, hDD22GF);
  // extrap(cctkGH, lambdaU0GF);
  // extrap(cctkGH, lambdaU1GF);
  // extrap(cctkGH, lambdaU2GF);
  // extrap(cctkGH, trKGF);
  // extrap(cctkGH, vetU0GF);
  // extrap(cctkGH, vetU1GF);
  // extrap(cctkGH, vetU2GF);
}

} // namespace NewRad
