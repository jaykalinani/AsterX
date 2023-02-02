/*! \file atmo.hxx
\brief Class definition representing artificial atmosphere.
*/

#ifndef ATMO_HXX
#define ATMO_HXX

#include "prims.hxx"
#include "cons.hxx"

namespace Con2PrimFactory {

/// Class representing an artificial atmosphere.
struct atmosphere {
  const CCTK_REAL rho_atmo;
  const CCTK_REAL eps_atmo;
  const CCTK_REAL ye_atmo;
  const CCTK_REAL press_atmo;
  const CCTK_REAL rho_cut;

  atmosphere(const atmosphere &) = default;
  atmosphere(CCTK_REAL rho_, CCTK_REAL eps_, CCTK_REAL Ye_, CCTK_REAL press_,
             CCTK_REAL rho_cut_);

  // Set prims to atmo
  CCTK_DEVICE CCTK_HOST void set(prim_vars &pv) const;

  // Set both prims and cons to atmo
  CCTK_DEVICE CCTK_HOST void set(prim_vars &pv, cons_vars &cv,
                                 const smat<CCTK_REAL, 3> &g) const;
};

} // namespace Con2PrimFactory
#endif
