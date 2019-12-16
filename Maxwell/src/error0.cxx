#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace Maxwell {

extern "C" void Maxwell_EstimateError(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_EstimateError;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> regrid_error_(cctkGH, regrid_error);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    const auto diffx = [&](const auto &var_, int i, int j, int k) {
      CCTK_REAL err = 0;
      for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a)
          err += fabs(var_(i + 1, j + a, k + b) - var_(i, j + a, k + b));
      return err;
    };
    const auto diffy = [&](const auto &var_, int i, int j, int k) {
      CCTK_REAL err = 0;
      for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a)
          err += fabs(var_(i + a, j + 1, k + b) - var_(i + a, j, k + b));
      return err;
    };
    const auto diffz = [&](const auto &var_, int i, int j, int k) {
      CCTK_REAL err = 0;
      for (int b = 0; b < 2; ++b)
        for (int a = 0; a < 2; ++a)
          err += fabs(var_(i + a, j + b, k + 1) - var_(i + a, j + b, k));
      return err;
    };

    const auto diff000 = [&](const auto &var_) {
      return diffx(var_, p.i, p.j, p.k) + diffy(var_, p.i, p.j, p.k) +
             diffz(var_, p.i, p.j, p.k);
    };
    const auto diff100 = [&](const auto &var_) {
      return diffx(var_, p.i - 1, p.j, p.k) + diffx(var_, p.i, p.j, p.k) +
             diffy(var_, p.i, p.j, p.k) + diffz(var_, p.i, p.j, p.k);
    };
    const auto diff010 = [&](const auto &var_) {
      return diffx(var_, p.i, p.j, p.k) + diffy(var_, p.i, p.j - 1, p.k) +
             diffy(var_, p.i, p.j, p.k) + diffz(var_, p.i, p.j, p.k);
    };
    const auto diff001 = [&](const auto &var_) {
      return diffx(var_, p.i, p.j, p.k) + diffy(var_, p.i, p.j, p.k) +
             diffz(var_, p.i, p.j, p.k - 1) + diffz(var_, p.i, p.j, p.k);
    };
    const auto diff011 = [&](const auto &var_) {
      return diffx(var_, p.i, p.j, p.k) + diffy(var_, p.i, p.j - 1, p.k) +
             diffy(var_, p.i, p.j, p.k) + diffz(var_, p.i, p.j, p.k - 1) +
             diffz(var_, p.i, p.j, p.k);
    };
    const auto diff101 = [&](const auto &var_) {
      return diffx(var_, p.i - 1, p.j, p.k) + diffx(var_, p.i, p.j, p.k) +
             diffy(var_, p.i, p.j, p.k) + diffz(var_, p.i, p.j, p.k - 1) +
             diffz(var_, p.i, p.j, p.k);
    };
    const auto diff110 = [&](const auto &var_) {
      return diffx(var_, p.i - 1, p.j, p.k) + diffx(var_, p.i, p.j, p.k) +
             diffy(var_, p.i, p.j - 1, p.k) + diffy(var_, p.i, p.j, p.k) +
             diffz(var_, p.i, p.j, p.k);
    };

    regrid_error_(p.I) = diff000(phi_) + diff100(ax_) + diff010(ay_) +
                         diff001(az_) + diff100(ex_) + diff010(ey_) +
                         diff001(ez_) + diff011(byz_) + diff101(bzx_) +
                         diff110(bxy_);
  });
}

} // namespace Maxwell
