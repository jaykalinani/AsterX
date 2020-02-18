#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_ConstraintBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ConstraintBoundaries;

  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCx_(cctkGH, ZtCx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCy_(cctkGH, ZtCy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCz_(cctkGH, ZtCz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_HC_(cctkGH, HC);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCx_(cctkGH, MtCx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCy_(cctkGH, MtCy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCz_(cctkGH, MtCz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_allC_(cctkGH, allC);

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_ZtCx_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_ZtCy_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_ZtCz_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_HC_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_MtCx_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_MtCy_(p.I) = 0; });
  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_MtCz_(p.I) = 0; });

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gf_allC_(p.I) = 0; });
}

} // namespace Z4c
