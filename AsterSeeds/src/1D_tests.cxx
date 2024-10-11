#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

#include "eos.hxx"
#include "eos_idealgas.hxx"
#include "seeds_utils.hxx"

namespace AsterSeeds {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void Tests1D_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_Tests1D_Initialize;
  DECLARE_CCTK_PARAMETERS;

  // For all the tests, the initial data EOS is ideal gas
  // Constructing the IG EOS object
  eos::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  const eos_idealgas eos_th(gl_gamma, particle_mass, rgeps, rgrho, rgye);
  const CCTK_REAL dummy_ye = 0.5;

  if (CCTK_EQUALS(test_case, "equilibrium")) {

    grid.loop_all<1, 1, 1>(grid.nghostzones,
                           [=] CCTK_HOST(const PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 rho(p.I) = 1.0;
                                 velx(p.I) = 0.0;
                                 vely(p.I) = 0.0;
                                 velz(p.I) = 0.0;
                                 press(p.I) = 1.0;
                                 eps(p.I) = eos_th.eps_from_valid_rho_press_ye(
                                     rho(p.I), press(p.I), dummy_ye);
                               });

    grid.loop_all<1, 0, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_x(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 1, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_y(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 0, 1>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_z(p.I) = 0.0;
                                                 });

  } else if (CCTK_EQUALS(test_case, "sound wave")) {
    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          vely(p.I) = 0.0;
          velz(p.I) = 0.0;
          press(p.I) = 1.0; // should add kinetic energy here
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all<1, 0, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_x(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 1, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_y(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 0, 1>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_z(p.I) = 0.0;
                                                 });

  } else if (CCTK_EQUALS(test_case, "Alfven wave")) {
    const CCTK_REAL A0 = 1.0;
    const CCTK_REAL va = 0.5;
    const CCTK_REAL k = 2 * M_PI;

    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          rho(p.I) = 1.0;
          velx(p.I) = 0.0;
          vely(p.I) = -va * A0 * cos(k * p.x);
          velz(p.I) = -va * A0 * sin(k * p.x);
          press(p.I) = 0.5; // should add kinetic energy here
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          Avec_x(p.I) = p.z * cos(k * p.x) - p.y * sin(k * p.x);
        });

    grid.loop_all<0, 1, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_y(p.I) = -p.z / 2.0;
                                                 });
    // CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_y(p.I) = p.x*sin(k*p.x) - p.z; });

    grid.loop_all<0, 0, 1>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_z(p.I) = p.y / 2.0;
                                                 });
    // CCTK_ATTRIBUTE_ALWAYS_INLINE { Avec_z(p.I) = p.y - p.x*cos(k*p.x); });

  } else if (CCTK_EQUALS(test_case, "shock tube")) {
    grid.loop_all<1, 1, 1>(grid.nghostzones,
                           [=] CCTK_HOST(const PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 if (p.x <= 0.0) {
                                   rho(p.I) = 2.0;
                                   velx(p.I) = 0.0;
                                   vely(p.I) = 0.0;
                                   velz(p.I) = 0.0;
                                   press(p.I) = 2.0;
                                 } else {
                                   rho(p.I) = 1.0;
                                   velx(p.I) = 0.0;
                                   vely(p.I) = 0.0;
                                   velz(p.I) = 0.0;
                                   press(p.I) = 1.0;
                                 }
                                 eps(p.I) = eos_th.eps_from_valid_rho_press_ye(
                                     rho(p.I), press(p.I), dummy_ye);
                               });

    grid.loop_all<1, 0, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_x(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 1, 0>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_y(p.I) = 0.0;
                                                 });

    grid.loop_all<0, 0, 1>(grid.nghostzones, [=] CCTK_HOST(const PointDesc &p)
                                                 CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                                   Avec_z(p.I) = 0.0;
                                                 });

  } else if (CCTK_EQUALS(test_case, "Balsara1") ||
             CCTK_EQUALS(test_case, "Balsara2") ||
             CCTK_EQUALS(test_case, "Balsara3") ||
             CCTK_EQUALS(test_case, "Balsara4") ||
             CCTK_EQUALS(test_case, "Balsara5")) {

    auto myNaN = numeric_limits<CCTK_REAL>::quiet_NaN();

    CCTK_REAL rhol = myNaN;
    CCTK_REAL vxl = myNaN;
    CCTK_REAL vyl = myNaN;
    CCTK_REAL vzl = myNaN;
    CCTK_REAL pressl = myNaN;
    CCTK_REAL Bxl = myNaN;
    CCTK_REAL Byl = myNaN;
    CCTK_REAL Bzl = myNaN;

    CCTK_REAL rhor = myNaN;
    CCTK_REAL vxr = myNaN;
    CCTK_REAL vyr = myNaN;
    CCTK_REAL vzr = myNaN;
    CCTK_REAL pressr = myNaN;
    CCTK_REAL Bxr = myNaN;
    CCTK_REAL Byr = myNaN;
    CCTK_REAL Bzr = myNaN;

    if (CCTK_EQUALS(test_case, "Balsara1")) {
      rhol = 1.0;
      vxl = 0.0;
      vyl = 0.0;
      vzl = 0.0;
      pressl = 1.0;
      Bxl = 0.5;
      Byl = 1.0;
      Bzl = 0.0;

      rhor = 0.125;
      vxr = 0.0;
      vyr = 0.0;
      vzr = 0.0;
      pressr = 0.1;
      Bxr = 0.5;
      Byr = -1.0;
      Bzr = 0.0;
    } else if (CCTK_EQUALS(test_case, "Balsara2")) {
      rhol = 1.0;
      vxl = 0.0;
      vyl = 0.0;
      vzl = 0.0;
      pressl = 30.0;
      Bxl = 5.0;
      Byl = 6.0;
      Bzl = 6.0;

      rhor = 1.0;
      vxr = 0.0;
      vyr = 0.0;
      vzr = 0.0;
      pressr = 1.0;
      Bxr = 5.0;
      Byr = 0.7;
      Bzr = 0.7;
    } else if (CCTK_EQUALS(test_case, "Balsara3")) {
      rhol = 1.0;
      vxl = 0.0;
      vyl = 0.0;
      vzl = 0.0;
      pressl = 1000.0;
      Bxl = 10.0;
      Byl = 7.0;
      Bzl = 7.0;

      rhor = 1.0;
      vxr = 0.0;
      vyr = 0.0;
      vzr = 0.0;
      pressr = 0.1;
      Bxr = 10.0;
      Byr = 0.7;
      Bzr = 0.7;
    } else if (CCTK_EQUALS(test_case, "Balsara4")) {
      rhol = 1.0;
      vxl = 0.999;
      vyl = 0.0;
      vzl = 0.0;
      pressl = 0.1;
      Bxl = 10.0;
      Byl = 7.0;
      Bzl = 7.0;

      rhor = 1.0;
      vxr = -0.999;
      vyr = 0.0;
      vzr = 0.0;
      pressr = 0.1;
      Bxr = 10.0;
      Byr = -7.0;
      Bzr = -7.0;
    } else if (CCTK_EQUALS(test_case, "Balsara5")) {
      rhol = 1.08;
      vxl = 0.4;
      vyl = 0.3;
      vzl = 0.2;
      pressl = 0.95;
      Bxl = 2.0;
      Byl = 0.3;
      Bzl = 0.3;

      rhor = 1.0;
      vxr = -0.45;
      vyr = -0.2;
      vzr = 0.2;
      pressr = 1.0;
      Bxr = 2.0;
      Byr = -0.7;
      Bzr = 0.5;
    } else {
      CCTK_ERROR("Balsara type not defined");
    }

    CCTK_REAL theta_x = myNaN;
    CCTK_REAL theta_y = myNaN;
    CCTK_REAL theta_z = myNaN;

    if (CCTK_EQUALS(shock_dir, "x")) {
      theta_x = 0.0;
      theta_y = 0.0;
      theta_z = 0.0;
    } else if (CCTK_EQUALS(shock_dir, "y")) {
      theta_x = 0.0;
      theta_y = 0.0;
      theta_z = 0.5 * M_PI;
    } else if (CCTK_EQUALS(shock_dir, "z")) {
      theta_x = 0.0;
      theta_y = 0.5 * M_PI;
      theta_z = 0.0;
    } else {
      theta_x = rotate_angle_x * M_PI;
      theta_y = rotate_angle_y * M_PI;
      theta_z = rotate_angle_z * M_PI;
    }

    // Rotation matrix R_{ij}
    const auto calc_R = [&](const CCTK_REAL th_x, const CCTK_REAL th_y,
                            const CCTK_REAL th_z) {
      const mat<CCTK_REAL, 3> R_x{
          1.0, 0.0, 0.0, 0.0, cos(th_x), -sin(th_x), 0.0, sin(th_x), cos(th_x)};
      const mat<CCTK_REAL, 3> R_y{cos(th_y),  0.0, sin(th_y), 0.0, 1.0, 0.0,
                                  -sin(th_y), 0.0, cos(th_y)};
      const mat<CCTK_REAL, 3> R_z{
          cos(th_z), -sin(th_z), 0.0, sin(th_z), cos(th_z), 0.0, 0.0, 0.0, 1.0};
      // R_{ij} = Rz_{ik} Ry_{kl} Rx_{lj}
      return mat<CCTK_REAL, 3>([&](int i, int j) ARITH_INLINE {
        return sum<3>([&](int k) ARITH_INLINE {
          return R_z(i, k) * sum<3>([&](int l) ARITH_INLINE {
                   return R_y(k, l) * R_x(l, j);
                 });
        });
      });
    };

    // Tranform vector
    const auto transform_vec = [&](const vec<CCTK_REAL, 3> u_vec,
                                   const mat<CCTK_REAL, 3> R_mat) {
      return vec<CCTK_REAL, 3>([&](int i) ARITH_INLINE {
        return sum<3>([&](int k)
                          ARITH_INLINE { return R_mat(i, k) * u_vec(k); });
      });
    };

    const auto R = calc_R(theta_x, theta_y, theta_z);
    const auto Rinv = calc_R(-theta_x, -theta_y, -theta_z);

    // B and v in current coordinates
    const auto oldBls = transform_vec(vec<CCTK_REAL, 3>{Bxl, Byl, Bzl}, Rinv);
    const auto oldBrs = transform_vec(vec<CCTK_REAL, 3>{Bxr, Byr, Bzr}, Rinv);
    const auto oldvls = transform_vec(vec<CCTK_REAL, 3>{vxl, vyl, vzl}, Rinv);
    const auto oldvrs = transform_vec(vec<CCTK_REAL, 3>{vxr, vyr, vzr}, Rinv);

    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto newxs = transform_vec(vec<CCTK_REAL, 3>{p.x, p.y, p.z}, R);

          if (newxs(0) <= 0.0) {
            rho(p.I) = rhol;
            velx(p.I) = oldvls(0);
            vely(p.I) = oldvls(1);
            velz(p.I) = oldvls(2);
            press(p.I) = pressl;

          } else {
            rho(p.I) = rhor;
            velx(p.I) = oldvrs(0);
            vely(p.I) = oldvrs(1);
            velz(p.I) = oldvrs(2);
            press(p.I) = pressr;
          }
          eps(p.I) = eos_th.eps_from_valid_rho_press_ye(rho(p.I), press(p.I),
                                                        dummy_ye);
        });

    grid.loop_all<1, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto newxs = transform_vec(vec<CCTK_REAL, 3>{p.x, p.y, p.z}, R);
          if (newxs(0) <= 0.0) {
            Avec_x(p.I) = oldBls(1) * (p.z);
          } else {
            Avec_x(p.I) = oldBrs(1) * (p.z);
          }
        });

    grid.loop_all<0, 1, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto newxs = transform_vec(vec<CCTK_REAL, 3>{p.x, p.y, p.z}, R);
          if (newxs(0) <= 0.0) {
            Avec_y(p.I) = oldBls(2) * (p.x);
          } else {
            Avec_y(p.I) = oldBrs(2) * (p.x);
          }
        });

    grid.loop_all<0, 0, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          const auto newxs = transform_vec(vec<CCTK_REAL, 3>{p.x, p.y, p.z}, R);
          if (newxs(0) <= 0.0) {
            Avec_z(p.I) = oldBls(0) * (p.y);
          } else {
            Avec_z(p.I) = oldBrs(0) * (p.y);
          }
        });

  } else {
    CCTK_ERROR("Test case not defined");
  }
}

} // namespace AsterSeeds
