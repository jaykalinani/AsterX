#include "physics.hxx"
#include "tensor.hxx"
#include "z4c_vars.hxx"

#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace Z4c {
using namespace Loop;
using namespace std;

extern "C" void Z4c_Constraints(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Z4c_Constraints;
  DECLARE_CCTK_PARAMETERS;

  for (int d = 0; d < 3; ++d)
    if (cctk_nghostzones[d] < 2)
      CCTK_VERROR("Need at least 2 ghost zones");

  const vec3<CCTK_REAL> dx{
      CCTK_DELTA_SPACE(0),
      CCTK_DELTA_SPACE(1),
      CCTK_DELTA_SPACE(2),
  };

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_chi_(cctkGH, chi);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxx_(cctkGH, gammatxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxy_(cctkGH, gammatxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatxz_(cctkGH, gammatxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatyy_(cctkGH, gammatyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatyz_(cctkGH, gammatyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_gammatzz_(cctkGH, gammatzz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Kh_(cctkGH, Kh);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxx_(cctkGH, Atxx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxy_(cctkGH, Atxy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atxz_(cctkGH, Atxz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atyy_(cctkGH, Atyy);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atyz_(cctkGH, Atyz);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Atzz_(cctkGH, Atzz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamtx_(cctkGH, Gamtx);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamty_(cctkGH, Gamty);
  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Gamtz_(cctkGH, Gamtz);

  const GF3D<const CCTK_REAL, 0, 0, 0> gf_Theta_(cctkGH, Theta);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCx_(cctkGH, ZtCx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCy_(cctkGH, ZtCy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_ZtCz_(cctkGH, ZtCz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_HC_(cctkGH, HC);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCx_(cctkGH, MtCx);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCy_(cctkGH, MtCy);
  const GF3D<CCTK_REAL, 0, 0, 0> gf_MtCz_(cctkGH, MtCz);

  const GF3D<CCTK_REAL, 0, 0, 0> gf_allC_(cctkGH, allC);

  loop_int<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // Load and calculate
    const z4c_vars<CCTK_REAL> vars(
        kappa1, kappa2, f_mu_L, f_mu_S, eta, //
        gf_chi_, gf_gammatxx_, gf_gammatxy_, gf_gammatxz_, gf_gammatyy_,
        gf_gammatyz_, gf_gammatzz_, gf_Kh_, gf_Atxx_, gf_Atxy_, gf_Atxz_,
        gf_Atyy_, gf_Atyz_, gf_Atzz_, gf_Gamtx_, gf_Gamty_, gf_Gamtz_,
        gf_Theta_,
        //
        // gf_alphaG_, gf_betaGx_, gf_betaGy_, gf_betaGz_,
        gf_chi_, gf_chi_, gf_chi_, gf_chi_,
        //
        p.I, dx);

    // Store
    vars.ZtC.store(gf_ZtCx_, gf_ZtCy_, gf_ZtCz_, p);
    gf_HC_(p.I) = vars.HC;
    vars.MtC.store(gf_MtCx_, gf_MtCy_, gf_MtCz_, p);
    gf_allC_(p.I) = vars.allC;
  });
}

} // namespace Z4c
