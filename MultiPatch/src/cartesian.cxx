#include "multipatch.hxx"

#include <cctk_Parameters.h>

namespace MultiPatch {

Patch makeCartesianPatch() {
  DECLARE_CCTK_PARAMETERS;
  const PatchConnection outer_boundary{true, -1};
  Patch patch0;
  patch0.ncells = {cartesian_ncells_i, cartesian_ncells_j, cartesian_ncells_k};
  patch0.xmin = {cartesian_xmin, cartesian_ymin, cartesian_zmin};
  patch0.xmax = {cartesian_xmax, cartesian_ymax, cartesian_zmax};
  patch0.is_cartesian = true;
  patch0.connections = {{outer_boundary, outer_boundary, outer_boundary},
                        {outer_boundary, outer_boundary, outer_boundary}};
  return patch0;
}

class CartesianPatchSystem final : PatchSystem {
public:
  CartesianPatchSystem() : PatchSystem({makeCartesianPatch()}) {}

  std::tuple<int, vec<CCTK_REAL, dim, UP> >
  global2local(const vec<CCTK_REAL, dim, UP> &x) const override {
    return std::make_tuple(0, x);
  }

  std::tuple<vec<CCTK_REAL, dim, UP>, vec<vec<CCTK_REAL, dim, DN>, dim, UP>,
             vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP> >
  d2local_dglobal2(int patch, const vec<CCTK_REAL, dim, UP> &a) const override {
    assert(patch >= 0 && patch < num_patches());
    return std::make_tuple(
        a, zero<vec<vec<CCTK_REAL, dim, DN>, dim, UP> >()(),
        zero<vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP> >()());
  }
};

} // namespace MultiPatch
