/*! \file c2p.hxx
\brief Defines a c2p
\author Jay Kalinani

c2p is effectively an interface to be used by different c2p implementations.

*/

#ifndef C2P_HXX
#define C2P_HXX

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <math.h>

#include "prims.hxx"
#include "cons.hxx"
#include "atmo.hxx"
#include "c2p_report.hxx"

#include "eos.hxx"
#include "eos_idealgas.hxx"

#include "c2p_utils.hxx"

namespace Con2PrimFactory {

using namespace AsterUtils;

constexpr CCTK_INT X = 0;
constexpr CCTK_INT Y = 1;
constexpr CCTK_INT Z = 2;

/* Abstract class c2p */
class c2p {
protected:
  /* The constructor must initialize the following variables */

  CCTK_INT maxIterations;
  CCTK_REAL tolerance;
  CCTK_REAL rho_strict;
  bool ye_lenient;
  CCTK_REAL v_lim;
  CCTK_REAL w_lim;
  CCTK_REAL vw_lim;
  CCTK_REAL Bsq_lim;
  atmosphere atmo;
  CCTK_REAL cons_error;
  bool use_zprim;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_Ssq_Exact(const vec<CCTK_REAL, 3> &mom,
                    const smat<CCTK_REAL, 3> &gup) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_Bsq_Exact(const vec<CCTK_REAL, 3> &B_up,
                    const smat<CCTK_REAL, 3> &glo) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
      get_BiSi_Exact(const vec<CCTK_REAL, 3> &Bvec, const vec<CCTK_REAL, 3> &mom) const;
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline vec<CCTK_REAL, 2>
      get_WLorentz_bsq_Seeds(const vec<CCTK_REAL, 3> &B_up, const vec<CCTK_REAL, 3> &v_up,
                             const smat<CCTK_REAL, 3> &glo) const;

  template <typename EOSType>
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  prims_floors_and_ceilings(const EOSType &eos_th, prim_vars &pv, const cons_vars &cv,
        const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const;
};

template <typename EOSType>
CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
c2p::prims_floors_and_ceilings(const EOSType &eos_th, prim_vars &pv, const cons_vars &cv,
                               const smat<CCTK_REAL, 3> &glo, c2p_report &rep) const {

  //const CCTK_REAL spatial_detg = calc_det(glo);
  //const CCTK_REAL sqrt_detg = sqrt(spatial_detg);

  // Lower velocity
  const vec<CCTK_REAL, 3> v_low = calc_contraction(glo, pv.vel);

  // ----------
  // Floor and ceiling for rho and velocity
  // Keeps pressure the same and changes eps
  // ----------

  // check if computed velocities are within the specified limit
  CCTK_REAL vsq_Sol = calc_contraction(v_low, pv.vel);
  CCTK_REAL sol_v = sqrt(vsq_Sol);
  if (sol_v > v_lim) {
    /*
    printf("(sol_v > v_lim) is true! \n");
    printf("sol_v, v_lim: %26.16e, %26.16e \n", sol_v, v_lim);
    */
    // add mass, keeps conserved density D
    pv.rho = cv.dens / w_lim;
    pv.eps = eos_th.eps_from_valid_rho_press_ye(pv.rho, pv.press, pv.Ye);
    pv.kappa =
        eos_th.kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    // if (pv.rho >= rho_strict) {
    //  rep.set_speed_limit({ sol_v, sol_v, sol_v });
    //  set_to_nan(pv, cv);
    //  return;
    //}
    pv.vel *= v_lim / sol_v;
    pv.w_lor = w_lim;
    // pv.eps = std::min(std::max(eos_th.rgeps.min, pv.eps),
    // eos_th.rgeps.max);
    // pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  if (pv.rho > rho_strict) {

    // remove mass, changes conserved density D
    pv.rho = rho_strict;
    pv.eps = eos_th.eps_from_valid_rho_press_ye(rho_strict, pv.press, pv.Ye);
    pv.kappa =
        eos_th.kappa_from_valid_rho_eps_ye(rho_strict, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  // ----------
  // Floor and ceiling for eps
  // Keeps rho the same and changes press
  // ----------

  // check the validity of the computed eps
  auto rgeps = eos_th.range_eps_from_valid_rho_ye(pv.rho, pv.Ye);
  if (pv.eps > rgeps.max) {
    // printf("(pv.eps > rgeps.max) is true, adjusting cons.. \n");
    // if (pv.rho >= rho_strict) {
    //  rep.set_range_eps(pv.eps); // sets adjust_cons to false by default
    //  rep.adjust_cons = true;
    //  set_to_nan(pv, cv);
    //  return;
    //}
    pv.eps = rgeps.max;
    pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.kappa = eos_th.kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  } else if (pv.eps < rgeps.min) {
    /*
    printf(
        "(pv.eps < rgeps.min) is true! pv.eps, rgeps.min: %26.16e, %26.16e
    \n",
        pv.eps, rgeps.min);
    printf(" Not adjusting cons.. \n");
    */
    // rep.set_range_eps(rgeps.min); // sets adjust_cons to true
    pv.eps = rgeps.min;
    pv.press = eos_th.press_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);
    pv.kappa = eos_th.kappa_from_valid_rho_eps_ye(pv.rho, pv.eps, pv.Ye);

    rep.adjust_cons = true;
  }

  // TODO: check validity for Ye

}

} // namespace Con2PrimFactory

#endif
