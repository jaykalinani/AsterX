#ifndef RECONX_RECONSTRUCT_HXX
#define RECONX_RECONSTRUCT_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include "monocentral.hxx"
#include "minmod.hxx"
#include "ppm.hxx"
// #include "eppm.hxx"
#include "wenoz.hxx"
#include "mp5.hxx"

namespace ReconX {

using namespace std;
using namespace Arith;

// enum class for different reconstruction routines

enum class reconstruction_t {
  Godunov,
  minmod,
  monocentral,
  ppm,
  // eppm,
  wenoz,
  mp5
};

inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE array<CCTK_REAL, 2>
reconstruct(const GF3D2<const CCTK_REAL> &gf_var, const PointDesc &p,
            const reconstruction_t &reconstruction, const int &dir,
            const bool &gf_is_rho, const GF3D2<const CCTK_REAL> &gf_press,
            const GF3D2<const CCTK_REAL> &gf_vel_dir,
            const reconstruct_params_t &reconstruct_params) {
  // Neighbouring "plus" and "minus" cell indices
  const auto Immm = p.I - 3 * p.DI[dir];
  const auto Imm = p.I - 2 * p.DI[dir];
  const auto Im = p.I - p.DI[dir];
  const auto Ip = p.I;
  const auto Ipp = p.I + p.DI[dir];
  const auto Ippp = p.I + 2 * p.DI[dir];

  switch (reconstruction) {

  case reconstruction_t::Godunov: {
    CCTK_REAL var_m = gf_var(Im);
    CCTK_REAL var_p = gf_var(Ip);
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::minmod: {
    // reconstructs values of Im and Ip at the common face between these
    // two cells
    CCTK_REAL var_slope_p = gf_var(Ipp) - gf_var(Ip);
    CCTK_REAL var_slope_c = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
    // reconstructed Im on its "plus/right" side
    CCTK_REAL var_m = gf_var(Im) + minmod(var_slope_c, var_slope_m) / 2;
    // reconstructed Ip on its "minus/left" side
    CCTK_REAL var_p = gf_var(Ip) - minmod(var_slope_p, var_slope_c) / 2;
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::monocentral: {
    // reconstructs values of Im and Ip at the common face between these
    // two cells
    // reconstructed Im on its "plus/right" side
    CCTK_REAL var_slope_p = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_slope_m = gf_var(Im) - gf_var(Imm);
    CCTK_REAL var_m = gf_var(Im) + monocentral(var_slope_p, var_slope_m) / 2;
    // reconstructed Ip on its "minus/left" side
    var_slope_p = gf_var(Ipp) - gf_var(Ip);
    var_slope_m = gf_var(Ip) - gf_var(Im);
    CCTK_REAL var_p = gf_var(Ip) - monocentral(var_slope_p, var_slope_m) / 2;
    return array<CCTK_REAL, 2>{var_m, var_p};
  }

  case reconstruction_t::ppm: {
    const array<const vect<int, dim>, 5> cells_Im = {Immm, Imm, Im, Ip, Ipp};
    const array<const vect<int, dim>, 5> cells_Ip = {Imm, Im, Ip, Ipp, Ippp};

    const array<CCTK_REAL, 2> rc_Im =
        ppm(gf_var, cells_Im, dir, gf_is_rho, gf_press, gf_vel_dir,
            reconstruct_params);
    const array<CCTK_REAL, 2> rc_Ip =
        ppm(gf_var, cells_Ip, dir, gf_is_rho, gf_press, gf_vel_dir,
            reconstruct_params);

    return array<CCTK_REAL, 2>{rc_Im.at(1), rc_Ip.at(0)};
  }

    // case reconstruction_t::eppm: {
    //   const array<const vect<int, dim>, 5> cells_Im = {Immm, Imm, Im, Ip,
    //   Ipp}; const array<const vect<int, dim>, 5> cells_Ip = {Imm, Im, Ip,
    //   Ipp, Ippp};

    //  const array<CCTK_REAL, 2> rc_Im =
    //      eppm(gf_var, cells_Im, dir, gf_is_rho, gf_press, gf_vel_dir,
    //           reconstruct_params);
    //  const array<CCTK_REAL, 2> rc_Ip =
    //      eppm(gf_var, cells_Ip, dir, gf_is_rho, gf_press, gf_vel_dir,
    //           reconstruct_params);

    //  return array<CCTK_REAL, 2>{rc_Im.at(1), rc_Ip.at(0)};
    //}

  case reconstruction_t::wenoz: {
    const array<const vect<int, dim>, 5> cells_Im = {Immm, Imm, Im, Ip, Ipp};
    const array<const vect<int, dim>, 5> cells_Ip = {Imm, Im, Ip, Ipp, Ippp};

    const array<CCTK_REAL, 2> rc_Im =
        wenoz(gf_var, cells_Im, reconstruct_params);
    const array<CCTK_REAL, 2> rc_Ip =
        wenoz(gf_var, cells_Ip, reconstruct_params);

    return array<CCTK_REAL, 2>{rc_Im.at(1), rc_Ip.at(0)};
  }

  case reconstruction_t::mp5: {
    // for the left cell, the plus side has sequence: Immm, Imm, Im, Ip, Ipp
    // for the left cell, the minus side has sequence: Ipp, Ip, Im, Im, Immm
    // here, we need the plus side
    const array<const vect<int, dim>, 5> cells_Im = {Immm, Imm, Im, Ip, Ipp};

    // for the right cell, the plus side has sequence: Imm, Im, Ip, Ipp, Ippp
    // for the right cell, the minus side has sequence: Ippp, Ipp, Ip, Im, Imm
    // here, we need the minus side

    const array<const vect<int, dim>, 5> cells_Ip = {Ippp, Ipp, Ip, Im, Imm};

    const CCTK_REAL rc_Im = mp5(gf_var, cells_Im, reconstruct_params);
    const CCTK_REAL rc_Ip = mp5(gf_var, cells_Ip, reconstruct_params);

    return array<CCTK_REAL, 2>{rc_Im, rc_Ip};
  }

  default:
    assert(0);
  }
}

} // namespace ReconX

#endif // RECONX_RECONSTRUCT_HXX
