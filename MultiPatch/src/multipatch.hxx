#ifndef MULTIPATCH_HXX
#define MULTIPATCH_HXX

#include <fixmath.hxx>
#include <cctk.h>

#include <mat.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <tuple>
#include <utility>
#include <vector>

namespace MultiPatch {
using namespace Arith;

constexpr int dim = 3;

struct PatchConnection {
  bool is_outer_boundary;
  int other_patch; // -1 if outer boundary
};

struct Patch {
  vect<int, dim> ncells;
  vect<CCTK_REAL, dim> xmin, xmax; // cell boundaries
  bool is_cartesian;               // Jacobian is trivial
  vect<vect<PatchConnection, dim>, 2> connections;
};

class PatchSystem {
public:
  std::vector<Patch> patches;
  int num_patches() const { return patches.size(); }

  PatchSystem() {}
  PatchSystem(std::vector<Patch> patches) : patches(std::move(patches)) {}

  virtual ~PatchSystem() {}

  virtual std::tuple<int, vec<CCTK_REAL, dim, UP> >
  global2local(const vec<CCTK_REAL, dim, UP> &x) const = 0;

  virtual vec<CCTK_REAL, dim, UP>
  local2global(int patch, const vec<CCTK_REAL, dim, UP> &a) const {
    const auto x_dx = dlocal_dglobal(patch, a);
    return std::get<0>(x_dx);
  }

  // Calculating global derivatives d/dx from local derivatives d/da
  // requires the Jacobian da/dx, and also its derivative d^2/dx^2 for
  // second derivatives

  // da/dx[i,j] = da[i] / dx[j]
  virtual std::tuple<vec<CCTK_REAL, dim, UP>,
                     vec<vec<CCTK_REAL, dim, DN>, dim, UP> >
  dlocal_dglobal(int patch, const vec<CCTK_REAL, dim, UP> &a) const {
    const auto x_dx_ddx = d2local_dglobal2(patch, a);
    return std::make_tuple(std::get<0>(x_dx_ddx), std::get<1>(x_dx_ddx));
  }

  // d^2a/dx^2[i,j,k] = d^2a[i] / dx^2[j,k]
  virtual std::tuple<vec<CCTK_REAL, dim, UP>,
                     vec<vec<CCTK_REAL, dim, DN>, dim, UP>,
                     vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP> >
  d2local_dglobal2(int patch, const vec<CCTK_REAL, dim, UP> &a) const = 0;
};

} // namespace MultiPatch

#endif // #ifndef MULTIPATCH_HXX
