#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace Maxwell {

extern "C" void Maxwell_Solve(CCTK_ARGUMENTS) { SolvePoisson(); }

extern "C" void Maxwell_UpdatePhi(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_UpdatePhi;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi1_(cctkGH, phi1);

  Loop::loop_int<0, 0, 0>(
      cctkGH, [&](const Loop::PointDesc &p) { phi_(p.I) -= phi1_(p.I); });
}
} // namespace Maxwell
