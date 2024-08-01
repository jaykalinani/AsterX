#ifndef EOS_3P_TABULATED3D_HXX
#define EOS_3P_TABULATED3D_HXX

#include <cctk.h>
#include <cmath>
#include <hdf5.h>
#include <mpi.h>
#include <string>
#include "../eos_3p.hxx"
#include "../utils/eos_brent.hxx"
#include "../utils/eos_linear_interp_ND.hxx"
#include "eos_readtable_scollapse.hxx"

namespace EOSX {

class eos_3p_tabulated3d : public eos_3p {

public:
  enum errors {
    NO_ERRORS = 0,
    RHO_TOO_HIGH,
    RHO_TOO_LOW,
    YE_TOO_HIGH,
    YE_TOO_LOW,
    TEMP_TOO_HIGH,
    TEMP_TOO_LOW,
    num_errors
  };

  enum EV {
    PRESS = 0,
    EPS,
    S,
    CS2,
    MUE,
    MUP,
    MUN,
    XA,
    XH,
    XN,
    XP,
    ABAR,
    ZBAR,
    NUM_VARS
  };

  CCTK_REAL gamma;
  range rgeps;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  init(const string &filename, range &rgeps_, const range &rgrho_, const range &rgye_) {
    eos_readtable_scollapse(filename);
    set_range_rho(rgrho_);
    set_range_ye(rgye_);
    rgeps_ = range_eps_from_valid_rho_ye(rgrho_.min, rgye_.min);
    // TODO: first compute temp as a function of rho, ye, and eps, and then
    // initialize its range For now, as dummy, we pass range of eps as range of
    // temp
    set_range_temp(rgeps_);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  logtemp_from_eps(const CCTK_REAL rho, CCTK_REAL &eps,
                   const CCTK_REAL ye) const {
    const auto lrho = log(rho);
    const auto leps = log(eps + energy_shift);

    // Get eps ranges
    const auto varsmin =
        interptable.interpolate<EV::EPS>(lrho, interptable.xmin<1>(), ye);
    const auto varsmax =
        interptable.interpolate<EV::EPS>(lrho, interptable.xmax<1>(), ye);

    if (leps <= varsmin[0]) {
      eps = exp(varsmin[0]) - energy_shift;
      return interptable.xmin<1>();
    }
    if (leps >= varsmax[0]) {
      eps = exp(varsmax[0]) - energy_shift;
      return interptable.xmax<1>();
    }

    // Root finding interface closure
    const auto func = [&](CCTK_REAL &lt) {
      const auto vars = interptable.interpolate<EV::EPS>(lrho, lt, ye);
      return leps - vars[0];
    };

    return zero_brent(interptable.xmin<1>(), interptable.xmax<1>(), 1.e-14,
                      func);
  }
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                               const CCTK_REAL ye) const {
    // error = checkbounds<true>(rho, temp, ye);
    const auto lrho = log(rho);
    const auto ltemp = log(temp);
    const auto vars = interptable.interpolate<EV::PRESS>(lrho, ltemp, ye);

    return exp(vars[0]);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                              const CCTK_REAL ye) const {
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = logtemp_from_eps(lrho, eps, ye);
    const auto vars = interptable.interpolate<EV::PRESS>(lrho, ltemp, ye);

    return exp(vars[0]);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                              const CCTK_REAL ye) const {

    assert(
        !"This routine should not be used. There is no monotonicity condition "
         "to enforce a succesfull inversion from eps(press). So you better "
         "rewrite your code to not require this call...");

    return 0;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                             const CCTK_REAL ye) const {

    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = log(temp);
    auto const vars = interptable.interpolate<EV::EPS>(lrho, ltemp, ye);
    const auto eps = exp(vars[0]) - energy_shift;
    return eps;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                              const CCTK_REAL ye) const {
    const auto lrho = log(rho);
    const auto ltemp = log(temp);
    const auto vars = interptable.interpolate<EV::CS2>(lrho, ltemp, ye);
    assert(vars[0] >= 0); // Soundspeed^2 should never ever be negative

    return sqrt(vars[0]);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                             const CCTK_REAL ye) const {
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = logtemp_from_eps(lrho, eps, ye);
    const auto vars = interptable.interpolate<EV::CS2>(lrho, ltemp, ye);
    assert(vars[0] >= 0); // Soundspeed^2 should never ever be negative

    return sqrt(vars[0]);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_valid_rho_eps_ye(const CCTK_REAL rho, CCTK_REAL &eps,
                             const CCTK_REAL ye) const {
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = logtemp_from_eps(lrho, eps, ye);

    return exp(ltemp);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_valid_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                                     CCTK_REAL &dpdeps, const CCTK_REAL rho,
                                     const CCTK_REAL eps,
                                     const CCTK_REAL ye) const {
    printf("press_derivs_from_valid_rho_eps_ye is not supported for now!");
    exit(EXIT_FAILURE);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                                 const CCTK_REAL ye) const {
    const CCTK_REAL lrho = log(rho);
    const CCTK_REAL ltemp = log(temp);
    const auto vars = interptable.interpolate<EV::S>(lrho, ltemp, ye);

    return vars[0];
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline range
  range_eps_from_valid_rho_ye(const CCTK_REAL rho, const CCTK_REAL ye) const {
    //    const CCTK_REAL lrho = log(rho);
    //    rgeps.min = interptable.interpolate<EV::EPS>(lrho,
    //    interptable.xmin<1>(), ye); rgeps.max =
    //    interptable.interpolate<EV::EPS>(lrho, interptable.xmax<1>(), ye);

    return rgeps;
  }
};
} // namespace EOSX

#endif
