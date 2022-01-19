#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace Maxwell {
using namespace Loop;

extern "C" void Maxwell_Boundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Boundaries;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < dim; ++d)
    assert(cctk_nghostzones[d] <= 1);

  const GF3D<CCTK_REAL, 0, 1, 1> dyz_(cctkGH, dyz);
  const GF3D<CCTK_REAL, 1, 0, 1> dzx_(cctkGH, dzx);
  const GF3D<CCTK_REAL, 1, 1, 0> dxy_(cctkGH, dxy);

  const GF3D<CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const GF3D<CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const GF3D<CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  loop_bnd<0, 1, 1>(cctkGH, [&](const PointDesc &p) {
    int f = p.NI[0] != 0 ? -1 : +1;
    dyz_(p.I) = f * dyz_(p.I - p.NI);
  });
  loop_bnd<1, 0, 1>(cctkGH, [&](const PointDesc &p) {
    int f = p.NI[1] != 0 ? -1 : +1;
    dzx_(p.I) = f * dzx_(p.I - p.NI);
  });
  loop_bnd<1, 1, 0>(cctkGH, [&](const PointDesc &p) {
    int f = p.NI[2] != 0 ? -1 : +1;
    dxy_(p.I) = f * dxy_(p.I - p.NI);
  });

  loop_bnd<0, 1, 1>(cctkGH, [&](const PointDesc &p) { byz_(p.I) = 0; });
  loop_bnd<1, 0, 1>(cctkGH, [&](const PointDesc &p) { bzx_(p.I) = 0; });
  loop_bnd<1, 1, 0>(cctkGH, [&](const PointDesc &p) { bxy_(p.I) = 0; });
}

} // namespace Maxwell
