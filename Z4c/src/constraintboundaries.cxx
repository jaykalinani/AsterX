#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

namespace Z4c {
using namespace Loop;

extern "C" void Z4c_ConstraintBoundaries(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_ConstraintBoundaries;

  const array<int, dim> indextype = {0, 0, 0};
  const array<int, dim> nghostzones = {cctk_nghostzones[0], cctk_nghostzones[1],
                                       cctk_nghostzones[2]};

  const GF3D1<CCTK_REAL> gf_ZtCx_(cctkGH, indextype, nghostzones, ZtCx);
  const GF3D1<CCTK_REAL> gf_ZtCy_(cctkGH, indextype, nghostzones, ZtCy);
  const GF3D1<CCTK_REAL> gf_ZtCz_(cctkGH, indextype, nghostzones, ZtCz);

  const GF3D1<CCTK_REAL> gf_HC_(cctkGH, indextype, nghostzones, HC);

  const GF3D1<CCTK_REAL> gf_MtCx_(cctkGH, indextype, nghostzones, MtCx);
  const GF3D1<CCTK_REAL> gf_MtCy_(cctkGH, indextype, nghostzones, MtCy);
  const GF3D1<CCTK_REAL> gf_MtCz_(cctkGH, indextype, nghostzones, MtCz);

  const GF3D1<CCTK_REAL> gf_allC_(cctkGH, indextype, nghostzones, allC);

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
