#include "multipatch.hxx"

#include <cctk_Parameters.h>

#include <cmath>

namespace MultiPatch {

class CubedSpherePatchSystem final : public PatchSystem {

  static constexpr int ncells_i = 10;
  static constexpr int ncells_j = 10;
  static constexpr int ncells_k = 10;

  static constexpr int ncells_rad = 10;

  static constexpr CCTK_REAL rmin = 1.0;
  static constexpr CCTK_REAL rmax = 3.0;

  static Patch makePatch0() {
    DECLARE_CCTK_PARAMETERS;
    Patch patch0;
    patch0.ncells = {ncells_i, ncells_j, ncells_k};
    patch0.xmin = {-rmin, -rmin, -rmin};
    patch0.xmax = {+rmax, +rmax, +rmax};
    patch0.is_cartesian = true;
    patch0.connections = {{{false, 1}, {false, 2}, {false, 3}},
                          {{false, 4}, {false, 5}, {false, 6}}};
    return patch0;
  }

  static Patch makePatch1() {
    DECLARE_CCTK_PARAMETERS;
    Patch patch1;
    patch1.ncells = {ncells_rad, ncells_j, ncells_k};
    patch1.xmin = {rmin, -1.0, -1.0};
    patch1.xmax = {rmax, +1.0, +1.0};
    patch1.is_cartesian = false;
    patch1.connections = {{{false, 0}, {false, 2}, {false, 3}},
                          {{true, -1}, {false, 5}, {false, 6}}};
    return patch1;
  }

public:
  CubedSpherePatchSystem() : PatchSystem({makePatch0(), makePatch1()}) {}

  std::tuple<int, vec<CCTK_REAL, dim, UP> >
  global2local(const vec<CCTK_REAL, dim, UP> &x) const override {
    if (abs(x(0)) <= rmin && abs(x(1)) <= rmin && abs(x(2)) <= rmin)
      return std::make_tuple(0, x);
    if (x(0) >= abs(x(1)) && x(0) >= abs(x(2)))
      return std::make_tuple(1, Arith::nan<vec<CCTK_REAL, dim, UP> >()());
    if (x(0) <= abs(x(1)) && x(0) <= abs(x(2)))
      return std::make_tuple(4, Arith::nan<vec<CCTK_REAL, dim, UP> >()());
    if (x(1) >= abs(x(0)) && x(1) >= abs(x(2)))
      return std::make_tuple(2, Arith::nan<vec<CCTK_REAL, dim, UP> >()());
    if (x(1) <= abs(x(0)) && x(1) <= abs(x(2)))
      return std::make_tuple(5, Arith::nan<vec<CCTK_REAL, dim, UP> >()());
    if (x(2) >= abs(x(0)) && x(2) >= abs(x(1)))
      return std::make_tuple(3, Arith::nan<vec<CCTK_REAL, dim, UP> >()());
    if (x(2) <= abs(x(0)) && x(2) <= abs(x(1)))
      return std::make_tuple(6, Arith::nan<vec<CCTK_REAL, dim, UP> >()());
    assert(0);
  }

  std::tuple<vec<CCTK_REAL, dim, UP>, vec<vec<CCTK_REAL, dim, DN>, dim, UP>,
             vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP> >
  d2local_dglobal2(int patch, const vec<CCTK_REAL, dim, UP> &a) const override {
    switch (patch) {
    case 0:
      return std::make_tuple(
          a, zero<vec<vec<CCTK_REAL, dim, DN>, dim, UP> >()(),
          zero<vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP> >()());
    case 1: {
      const auto rho = a(0);
      const auto mu_y = a(1);
      const auto mu_z = a(2);

      // Conditions:
      // - \mu \in [-1; +1]
      // - for rho = rmin, \mu lie on a cube at rmin
      // - for rho = rmax, \mu lie on a sphere at rmax

      // for rho = rmin:
      const auto x0 = rmin;
      const auto y0 = rmin * mu_y;
      const auto z0 = rmin * mu_z;
      const auto r0 = rmin * sqrt(pow2(x0) + pow2(y0) + pow2(z0));
      // y = r \cos \theta
      // z = r \cos \theta
      const auto theta_y = acos(y0 / r0);
      const auto theta_z = acos(z0 / r0);
      // for rho = rmax:
      const auto r1 = rmax;
      const auto r = lincom(rmin, r0, rmax, r1, rho);

      const auto y = r * cos(theta_y);
      const auto z = r * cos(theta_z);
      const auto x = sqrt(pow2(r) - pow2(y) - pow2(z));

      return std::make_tuple(
          vec<CCTK_REAL, dim, UP>{x, y, z},
          Arith::nan<vec<vec<CCTK_REAL, dim, DN>, dim, UP> >()(),
          Arith::nan<vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP> >()());
    }
    // case 2:
    //   ...;
    // case 3:
    //   ...;
    // case 4:
    //   ...;
    // case 5:
    //   ...;
    // case 6:
    //   ...;
    default:
      assert(0);
    }
  }
};

std::unique_ptr<PatchSystem> SetupCubedSphere() {
  return std::make_unique<CubedSpherePatchSystem>();
}

} // namespace MultiPatch
