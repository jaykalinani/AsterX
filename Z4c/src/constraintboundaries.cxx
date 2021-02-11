#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_ConstraintBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ConstraintBoundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const GF3D2layout layout(cctkGH, indextype);

  const GF3D2<CCTK_REAL> gf_ZtCx_(&layout, ZtCx);
  const GF3D2<CCTK_REAL> gf_ZtCy_(&layout, ZtCy);
  const GF3D2<CCTK_REAL> gf_ZtCz_(&layout, ZtCz);

  const GF3D2<CCTK_REAL> gf_HC_(&layout, HC);

  const GF3D2<CCTK_REAL> gf_MtCx_(&layout, MtCx);
  const GF3D2<CCTK_REAL> gf_MtCy_(&layout, MtCy);
  const GF3D2<CCTK_REAL> gf_MtCz_(&layout, MtCz);

  const GF3D2<CCTK_REAL> gf_allC_(&layout, allC);

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
