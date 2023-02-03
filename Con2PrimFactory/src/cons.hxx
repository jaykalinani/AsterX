/*! \file cons.hxx
\brief Class definitions representing conserved mhd variables.
*/

#ifndef CONS_HXX
#define CONS_HXX

#include "prims.hxx"

namespace Con2PrimFactory {

// Structure to represent conserved variables for mhd
struct cons_vars {
  CCTK_REAL dens;
  vec<CCTK_REAL, 3> mom;
  CCTK_REAL tau;
  CCTK_REAL dYe;
  vec<CCTK_REAL, 3> dBvec;

  // Default constructor, no initialization.
  cons_vars() = default;

  // Construct from single variables.
  CCTK_HOST CCTK_DEVICE cons_vars(CCTK_REAL dens_, vec<CCTK_REAL, 3> mom_, CCTK_REAL tau_,
            CCTK_REAL dYe_, vec<CCTK_REAL, 3> dBvec_)
      : dens{dens_}, mom{mom_}, tau{tau_}, dYe{dYe_}, dBvec{dBvec_} {}

  // Compute conserved variables from primitives and 3-metric
  CCTK_HOST CCTK_DEVICE void from_prim(const smat<CCTK_REAL, 3> &g,
                                       const CCTK_REAL &lapse,
                                       const vec<CCTK_REAL, 3> &beta_up,
                                       const prim_vars &pv);

  // Copy all members into single variables
  CCTK_HOST CCTK_DEVICE void scatter(CCTK_REAL &dens_, CCTK_REAL &momx_,
                                     CCTK_REAL &momy_, CCTK_REAL &momz_,
                                     CCTK_REAL &tau_, CCTK_REAL &dYe_,
                                     CCTK_REAL &dBvecx_, CCTK_REAL &dBvecy_,
                                     CCTK_REAL &dBvecz_) const;

  // Set all data to NAN
  void set_to_nan();
};

} // namespace Con2PrimFactory
#endif
