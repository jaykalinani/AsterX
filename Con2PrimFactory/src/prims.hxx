/*! \file prims.hxx
\brief Class definitions representing primitive mhd variables.
*/

#ifndef PRIMS_HXX
#define PRIMS_HXX

#include "c2p_utils.hxx"

namespace Con2PrimFactory {

struct prim_vars {
  CCTK_REAL rho;
  CCTK_REAL eps;
  CCTK_REAL Ye;
  CCTK_REAL press;
  vec<CCTK_REAL, 3> vel;
  CCTK_REAL w_lor;
  vec<CCTK_REAL, 3> Bvec;
  vec<CCTK_REAL, 3> E;

  /// Default constructor. Leaves all members uninitialized.
  CCTK_HOST CCTK_DEVICE prim_vars() = default;

  /// Construct from single variables.
  CCTK_HOST CCTK_DEVICE prim_vars(CCTK_REAL rho_, CCTK_REAL eps_, CCTK_REAL ye_,
                                  CCTK_REAL press_, vec<CCTK_REAL, 3> vel_,
                                  CCTK_REAL w_lor_, vec<CCTK_REAL, 3> Bvec_)
      : rho(rho_), eps(eps_), Ye(ye_), press(press_), vel(vel_), w_lor(w_lor_),
        Bvec(Bvec_){};

  /// Copy the members into single variables.
  CCTK_HOST CCTK_DEVICE void scatter(CCTK_REAL &rho_, CCTK_REAL &eps_,
                                     CCTK_REAL &ye_, CCTK_REAL &press_,
                                     CCTK_REAL &velx_, CCTK_REAL &vely_,
                                     CCTK_REAL &velz_, CCTK_REAL &w_lor_,
                                     CCTK_REAL &Bvec_x_, CCTK_REAL &Bvec_y_,
                                     CCTK_REAL &Bvec_z_, CCTK_REAL &E_x_,
                                     CCTK_REAL &E_y_, CCTK_REAL &E_z_) const;

  /// Set all data to NAN
  CCTK_HOST CCTK_DEVICE void set_to_nan();
};

} // namespace Con2PrimFactory
#endif
