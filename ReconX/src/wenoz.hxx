#ifndef RECONX_WENOZ_HXX
#define RECONX_WENOZ_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>

#include "reconx_utils.hxx"

namespace ReconX {

using namespace std;
using namespace Arith;
using namespace Loop;

/* WENO-Z: Performs the reconstruction of a given variable using 5th order
 * WENO-Z method based on Borges et al., ''An improved weighted essentially
 * non-oscillatory scheme for hyperbolic conservation laws'', 2008) */
// Also, see the Spritz code for details

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST array<CCTK_REAL, 2>
wenoz(const GF3D2<const CCTK_REAL> &gf_var,
      const array<const vect<int, dim>, 5> &cells,
      const reconstruct_params_t &reconstruct_params) {

  // Unpack all WENO-Z parameters
  const CCTK_REAL weno_eps = reconstruct_params.weno_eps;

  // Unpack all cells in the stencil
  const auto &Imm = cells[0];
  const auto &Im = cells[1];
  const auto &I = cells[2];
  const auto &Ip = cells[3];
  const auto &Ipp = cells[4];

  // Grid function at neighboring cells
  const CCTK_REAL &gf_Imm = gf_var(Imm);
  const CCTK_REAL &gf_Im = gf_var(Im);
  const CCTK_REAL &gf_I = gf_var(I);
  const CCTK_REAL &gf_Ip = gf_var(Ip);
  const CCTK_REAL &gf_Ipp = gf_var(Ipp);

  // Computing the smoothness indicators (Borges et al. 2008)
  const vec<CCTK_REAL, 3> betaZ{
      (13.0 / 12.0) * pow2(gf_Imm - 2.0 * gf_Im + gf_I) +
          (1.0 / 4.0) * pow2(gf_Imm - 4.0 * gf_Im + 3.0 * gf_I),

      (13.0 / 12.0) * pow2(gf_Im - 2.0 * gf_I + gf_Ip) +
          (1.0 / 4.0) * pow2(gf_Im - gf_Ip),

      (13.0 / 12.0) * pow2(gf_I - 2.0 * gf_Ip + gf_Ipp) +
          (1.0 / 4.0) * pow2(3.0 * gf_I - 4.0 * gf_Ip + gf_Ipp)};

  // Defining tau5 based on eq. 25 of (Borges et al. 2008)
  const CCTK_REAL tau5 = abs(betaZ(0) - betaZ(2));

  // Unnormalized weights
  // Optimal weights are chosen according to Del Zanna et al. 2007

  vec<CCTK_REAL, 3> aux_alphaZ;
  aux_alphaZ(0) = 1.0 + tau5 / (betaZ(0) + weno_eps);
  aux_alphaZ(1) = 1.0 + tau5 / (betaZ(1) + weno_eps);
  aux_alphaZ(2) = 1.0 + tau5 / (betaZ(2) + weno_eps);
  // const vec<CCTK_REAL, 3> wt = {5.0 / 16.0, 10.0 / 16.0, 1.0 / 16.0};
  // Original weights as suggested in (Borges et al. 2008)
  const vec<CCTK_REAL, 3> wt = {3.0 / 10.0, 3.0 / 5.0, 1.0 / 10.0};

  vec<vec<CCTK_REAL, 2>, 3> alphaZ;

  // for minus side
  alphaZ(0)(0) = wt(0) * aux_alphaZ(0);
  alphaZ(1)(0) = wt(1) * aux_alphaZ(1);
  alphaZ(2)(0) = wt(2) * aux_alphaZ(2);

  // for plus side
  alphaZ(0)(1) = wt(2) * aux_alphaZ(0);
  alphaZ(1)(1) = wt(1) * aux_alphaZ(1);
  alphaZ(2)(1) = wt(0) * aux_alphaZ(2);

  // Normalized weights for the reconstruction (Del Zanna et al. 2007)
  const vec<CCTK_REAL, 2> omega_denom([&](int j) ARITH_INLINE {
    return alphaZ(0)(j) + alphaZ(1)(j) + alphaZ(2)(j);
  });

  const vec<vec<CCTK_REAL, 2>, 3> omegaZ([&](int j) ARITH_INLINE {
    return vec<CCTK_REAL, 2>(
        [&](int f) ARITH_INLINE { return alphaZ(j)(f) / omega_denom(f); });
  });

  // Reconstruct cell-centered variable to left (minus) and right (plus) cell
  // interfaces

  /* Spritz Weights:
  const CCTK_REAL var_m =
      (omegaZ(2)(0) / 8.0) * (3.0 * gf_Ipp - 10.0 * gf_Ip + 15.0 * gf_I) +
      (omegaZ(1)(0) / 8.0) * (-1.0 * gf_Ip + 6.0 * gf_I + 3.0 * gf_Im) +
      (omegaZ(0)(0) / 8.0) * (3.0 * gf_I + 6.0 * gf_Im - 1.0 * gf_Imm);

  const CCTK_REAL var_p =
      (omegaZ(0)(1) / 8.0) * (3.0 * gf_Imm - 10.0 * gf_Im + 15.0 * gf_I) +
      (omegaZ(1)(1) / 8.0) * (-1.0 * gf_Im + 6.0 * gf_I + 3.0 * gf_Ip) +
      (omegaZ(2)(1) / 8.0) * (3.0 * gf_I + 6.0 * gf_Ip - 1.0 * gf_Ipp);
  */

  // GRHydro Weights:
  const CCTK_REAL var_m =
      (omegaZ(2)(0) / 6.0) * (2.0 * gf_Ipp - 7.0 * gf_Ip + 11.0 * gf_I) +
      (omegaZ(1)(0) / 6.0) * (-1.0 * gf_Ip + 5.0 * gf_I + 2.0 * gf_Im) +
      (omegaZ(0)(0) / 6.0) * (2.0 * gf_I + 5.0 * gf_Im - 1.0 * gf_Imm);

  const CCTK_REAL var_p =
      (omegaZ(0)(1) / 6.0) * (2.0 * gf_Imm - 7.0 * gf_Im + 11.0 * gf_I) +
      (omegaZ(1)(1) / 6.0) * (-1.0 * gf_Im + 5.0 * gf_I + 2.0 * gf_Ip) +
      (omegaZ(2)(1) / 6.0) * (2.0 * gf_I + 5.0 * gf_Ip - 1.0 * gf_Ipp);

  const array<CCTK_REAL, 2> var_rc = {var_m, var_p};
  return var_rc;
}

} // namespace ReconX

#endif // RECONX_WENOZ_HXX
