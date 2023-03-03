#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <Vc/Vc>

#include <cassert>
#include <cmath>
#include <iostream>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

extern "C" void HydroToyCarpetX_Evolve(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Evolve;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_REAL dt = CCTK_DELTA_TIME;
  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const Loop::GF3D1<const CCTK_REAL> rho_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                            rho_p);
  const Loop::GF3D1<const CCTK_REAL> momx_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momx_p);
  const Loop::GF3D1<const CCTK_REAL> momy_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momy_p);
  const Loop::GF3D1<const CCTK_REAL> momz_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             momz_p);
  const Loop::GF3D1<const CCTK_REAL> etot_p_(cctkGH, {1, 1, 1}, {1, 1, 1},
                                             etot_p);

  const Loop::GF3D1<const CCTK_REAL> fxrho_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                            fxrho);
  const Loop::GF3D1<const CCTK_REAL> fxmomx_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                             fxmomx);
  const Loop::GF3D1<const CCTK_REAL> fxmomy_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                             fxmomy);
  const Loop::GF3D1<const CCTK_REAL> fxmomz_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                             fxmomz);
  const Loop::GF3D1<const CCTK_REAL> fxetot_(cctkGH, {0, 1, 1}, {0, 0, 0},
                                             fxetot);

  const Loop::GF3D1<const CCTK_REAL> fyrho_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                            fyrho);
  const Loop::GF3D1<const CCTK_REAL> fymomx_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                             fymomx);
  const Loop::GF3D1<const CCTK_REAL> fymomy_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                             fymomy);
  const Loop::GF3D1<const CCTK_REAL> fymomz_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                             fymomz);
  const Loop::GF3D1<const CCTK_REAL> fyetot_(cctkGH, {1, 0, 1}, {0, 0, 0},
                                             fyetot);

  const Loop::GF3D1<const CCTK_REAL> fzrho_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                            fzrho);
  const Loop::GF3D1<const CCTK_REAL> fzmomx_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                             fzmomx);
  const Loop::GF3D1<const CCTK_REAL> fzmomy_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                             fzmomy);
  const Loop::GF3D1<const CCTK_REAL> fzmomz_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                             fzmomz);
  const Loop::GF3D1<const CCTK_REAL> fzetot_(cctkGH, {1, 1, 0}, {0, 0, 0},
                                             fzetot);

  const Loop::GF3D1<CCTK_REAL> rho_(cctkGH, {1, 1, 1}, {1, 1, 1}, rho);
  const Loop::GF3D1<CCTK_REAL> momx_(cctkGH, {1, 1, 1}, {1, 1, 1}, momx);
  const Loop::GF3D1<CCTK_REAL> momy_(cctkGH, {1, 1, 1}, {1, 1, 1}, momy);
  const Loop::GF3D1<CCTK_REAL> momz_(cctkGH, {1, 1, 1}, {1, 1, 1}, momz);
  const Loop::GF3D1<CCTK_REAL> etot_(cctkGH, {1, 1, 1}, {1, 1, 1}, etot);

  // Transport
  // dt rho + d_i (rho vel^i) = 0
  // dt mom_j + d_i (mom_j vel^i) = 0
  // dt etot + d_i (etot vel^i) = 0

  const CCTK_REAL dt_dx = dt / dx;
  const CCTK_REAL dt_dy = dt / dy;
  const CCTK_REAL dt_dz = dt / dz;

  if (false) {
    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      auto calcupdate{[&](auto &fx, auto &fy, auto &fz) {
        return dt_dx * (fx(p.I + p.DI(0)) - fx(p.I)) +
               dt_dy * (fy(p.I + p.DI(1)) - fy(p.I)) +
               dt_dz * (fz(p.I + p.DI(2)) - fz(p.I));
      }};

      rho_(p.I) = rho_p_(p.I) - calcupdate(fxrho_, fyrho_, fzrho_);

      // momx_(p.I) = momx_p_(p.I) - calcupdate(fxmomx_, fymomx_, fzmomx_);
      // momy_(p.I) = momy_p_(p.I) - calcupdate(fxmomy_, fymomy_, fzmomy_);
      // momz_(p.I) = momz_p_(p.I) - calcupdate(fxmomz_, fymomz_, fzmomz_);

      // etot_(p.I) = etot_p_(p.I) - calcupdate(fxetot_, fyetot_, fzetot_);
    });
  }

  if (false) {
    const auto calcupdate{[&](auto I, auto &fx, auto &fy, auto &fz) {
      const Loop::vect<int, dim> DI{{1, 0, 0}};
      const Loop::vect<int, dim> DJ{{0, 1, 0}};
      const Loop::vect<int, dim> DK{{0, 0, 1}};
      return dt_dx * (fx(I + DI) - fx(I)) + dt_dy * (fy(I + DJ) - fy(I)) +
             dt_dz * (fz(I + DK) - fz(I));
    }};

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      rho_(p.I) = rho_p_(p.I) - calcupdate(p.I, fxrho_, fyrho_, fzrho_);
    });

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      momx_(p.I) = momx_p_(p.I) - calcupdate(p.I, fxmomx_, fymomx_, fzmomx_);
    });
    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      momy_(p.I) = momy_p_(p.I) - calcupdate(p.I, fxmomy_, fymomy_, fzmomy_);
    });
    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      momz_(p.I) = momz_p_(p.I) - calcupdate(p.I, fxmomz_, fymomz_, fzmomz_);
    });

    Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
      etot_(p.I) = etot_p_(p.I) - calcupdate(p.I, fxetot_, fyetot_, fzetot_);
    });
  }

  if (false) {
    Loop::vect<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = 1;
      imax[d] = cctk_lsh[d] - 1;
    }

    auto calcupdate{[&](auto I, auto &fx, auto &fy, auto &fz) {
      const Loop::vect<int, dim> DI{{1, 0, 0}};
      const Loop::vect<int, dim> DJ{{0, 1, 0}};
      const Loop::vect<int, dim> DK{{0, 0, 1}};
      return dt_dx * (fx(I + DI) - fx(I)) + dt_dy * (fy(I + DJ) - fy(I)) +
             dt_dz * (fz(I + DK) - fz(I));
    }};

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; ++i) {
          const Loop::vect<int, dim> I{{i, j, k}};
          rho_(I) = rho_p_(I) - calcupdate(I, fxrho_, fyrho_, fzrho_);
        }
      }
    }

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; ++i) {
          const Loop::vect<int, dim> I{{i, j, k}};
          momx_(I) = momx_p_(I) - calcupdate(I, fxmomx_, fymomx_, fzmomx_);
        }
      }
    }
    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; ++i) {
          const Loop::vect<int, dim> I{{i, j, k}};
          momy_(I) = momy_p_(I) - calcupdate(I, fxmomy_, fymomy_, fzmomy_);
        }
      }
    }
    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; ++i) {
          const Loop::vect<int, dim> I{{i, j, k}};
          momz_(I) = momz_p_(I) - calcupdate(I, fxmomz_, fymomz_, fzmomz_);
        }
      }
    }

    for (int k = imin[2]; k < imax[2]; ++k) {
      for (int j = imin[1]; j < imax[1]; ++j) {
        for (int i = imin[0]; i < imax[0]; ++i) {
          const Loop::vect<int, dim> I{{i, j, k}};
          etot_(I) = etot_p_(I) - calcupdate(I, fxetot_, fyetot_, fzetot_);
        }
      }
    }
  }

  if (false) {
    Loop::vect<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = 1;
      imax[d] = cctk_lsh[d] - 1;
    }

    const auto calcupdate{
        [&](auto &u, auto &u_p, auto &fx, auto &fy, auto &fz) {
          const Loop::vect<int, dim> DI{{1, 0, 0}};
          const Loop::vect<int, dim> DJ{{0, 1, 0}};
          const Loop::vect<int, dim> DK{{0, 0, 1}};
          for (int k = imin[2]; k < imax[2]; ++k) {
            for (int j = imin[1]; j < imax[1]; ++j) {
              for (int i = imin[0]; i < imax[0]; ++i) {
                const Loop::vect<int, dim> I{{i, j, k}};
                u(I) = u_p(I) - (dt_dx * (fx(I + DI) - fx(I)) +
                                 dt_dy * (fy(I + DJ) - fy(I)) +
                                 dt_dz * (fz(I + DK) - fz(I)));
              }
            }
          }
        }};

    calcupdate(rho_, rho_p_, fxrho_, fyrho_, fzrho_);
    calcupdate(momx_, momx_p_, fxmomx_, fymomx_, fzmomx_);
    calcupdate(momy_, momy_p_, fxmomy_, fymomy_, fzmomy_);
    calcupdate(momz_, momz_p_, fxmomz_, fymomz_, fzmomz_);
    calcupdate(etot_, etot_p_, fxetot_, fyetot_, fzetot_);
  }

  if (true) {

    Loop::vect<int, dim> ash;
    for (int d = 0; d < dim; ++d)
      ash[d] = cctk_ash[d];

    Loop::vect<int, dim> str[8];
    for (int itype = 0b000; itype <= 0b111; ++itype) {
      str[itype][0] = 1;
      str[itype][1] = str[itype][0] * (ash[0] + ((itype & 0b100) == 0));
      str[itype][2] = str[itype][1] * (ash[1] + ((itype & 0b010) == 0));
    }

    Loop::vect<int, dim> imin, imax;
    for (int d = 0; d < dim; ++d) {
      imin[d] = 1;
      imax[d] = cctk_lsh[d] - 1;
    }

    typedef Vc::native_simd<CCTK_REAL> CCTK_REAL_VEC;
    constexpr int VS = CCTK_REAL_VEC::size();
    assert((imax[0] - imin[0]) % VS == 0);

    auto vec_dt_dx = CCTK_REAL_VEC(dt_dx);
    auto vec_dt_dy = CCTK_REAL_VEC(dt_dy);
    auto vec_dt_dz = CCTK_REAL_VEC(dt_dz);

    const auto calcupdate{[&](auto &u_, auto &u_p_, auto &fx_, auto &fy_,
                              auto &fz_) {
      for (int k = imin[2]; k < imax[2]; ++k) {
        for (int j = imin[1]; j < imax[1]; ++j) {
          int idx111_i0 = str[0b111][1] * j + str[0b111][2] * k;
          int idx011_i0 = str[0b011][1] * j + str[0b011][2] * k;
          int idx101_i0 = str[0b101][1] * j + str[0b101][2] * k;
          int idx110_i0 = str[0b110][1] * j + str[0b110][2] * k;
          for (int i = imin[0]; i < imax[0]; i += VS) {
            int idx111 = idx111_i0 + i;
            int idx011 = idx011_i0 + i;
            int idx101 = idx101_i0 + i;
            int idx110 = idx110_i0 + i;
            CCTK_REAL_VEC u_p(&u_p_[idx111], Vc::element_aligned);
            CCTK_REAL_VEC fx1(&fx_[idx011 + str[0b011][0]],
                              Vc::element_aligned);
            CCTK_REAL_VEC fx0(&fx_[idx011], Vc::element_aligned);
            CCTK_REAL_VEC fy1(&fy_[idx101 + str[0b101][1]],
                              Vc::element_aligned);
            CCTK_REAL_VEC fy0(&fy_[idx101], Vc::element_aligned);
            CCTK_REAL_VEC fz1(&fz_[idx110 + str[0b110][2]],
                              Vc::element_aligned);
            CCTK_REAL_VEC fz0(&fz_[idx110], Vc::element_aligned);
            auto u = u_p - (vec_dt_dx * (fx1 - fx0) + vec_dt_dy * (fy1 - fy0) +
                            vec_dt_dz * (fz1 - fz0));
            u.copy_to(&u_[idx111], Vc::element_aligned);
          }
        }
      }
    }};

    calcupdate(rho, rho_p, fxrho, fyrho, fzrho);
    calcupdate(momx, momx_p, fxmomx, fymomx, fzmomx);
    calcupdate(momy, momy_p, fxmomy, fymomy, fzmomy);
    calcupdate(momz, momz_p, fxmomz, fymomz, fzmomz);
    calcupdate(etot, etot_p, fxetot, fyetot, fzetot);
  }
}

} // namespace HydroToyCarpetX
