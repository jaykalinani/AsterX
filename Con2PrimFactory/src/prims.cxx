/*! \file prims.cxx
\brief Implementation of classes for primitive mhd variables.
*/

#include "prims.hxx"
#include <limits>

namespace Con2PrimFactory {

CCTK_HOST CCTK_DEVICE void
prim_vars::scatter(CCTK_REAL &rho_, CCTK_REAL &eps_, CCTK_REAL &ye_,
                   CCTK_REAL &press_, CCTK_REAL &velx_, CCTK_REAL &vely_,
                   CCTK_REAL &velz_, CCTK_REAL &w_lor_, CCTK_REAL &Bvec_x_,
                   CCTK_REAL &Bvec_y_, CCTK_REAL &Bvec_z_, CCTK_REAL &E_x_,
                   CCTK_REAL &E_y_, CCTK_REAL &E_z_) const {
  rho_ = rho;
  eps_ = eps;
  ye_ = Ye;
  press_ = press;
  velx_ = vel(0);
  vely_ = vel(1);
  velz_ = vel(2);
  w_lor_ = w_lor;
  Bvec_x_ = Bvec(0);
  Bvec_y_ = Bvec(1);
  Bvec_z_ = Bvec(2);
  E_x_ = E(0);
  E_y_ = E(1);
  E_z_ = E(2);
}

CCTK_HOST CCTK_DEVICE void prim_vars::set_to_nan() {
  rho = eps = Ye = press = vel(0) = vel(1) = vel(2) = w_lor = Bvec(0) =
      Bvec(1) = Bvec(2) = E(0) = E(1) = E(2) =
          std::numeric_limits<CCTK_REAL>::quiet_NaN();
}

} // namespace Con2PrimFactory
