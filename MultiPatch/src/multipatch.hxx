#ifndef MULTIPATCH_HXX
#define MULTIPATCH_HXX

#include <fixmath.hxx>
#include <cctk.h>

#include <loop.hxx>
#include <mat.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace MultiPatch {
using namespace Arith;

constexpr int dim = 3;

struct PatchFace {
  bool is_outer_boundary;
  int other_patch; // -1 if outer boundary
};

struct Patch {
  vect<int, dim> ncells;
  vect<CCTK_REAL, dim> xmin, xmax; // cell boundaries
  bool is_cartesian;               // Jacobian is trivial
  vect<vect<PatchFace, dim>, 2> faces;
};

struct PatchTransformations {
  // Cartesian
  const CCTK_REAL cartesian_xmax;
  const CCTK_REAL cartesian_xmin;
  const CCTK_REAL cartesian_ymax;
  const CCTK_REAL cartesian_ymin;
  const CCTK_REAL cartesian_zmax;
  const CCTK_REAL cartesian_zmin;
  const int cartesian_ncells_i;
  const int cartesian_ncells_j;
  const int cartesian_ncells_k;

  // Cubed sphere
  const CCTK_REAL cubed_sphere_rmin;
  const CCTK_REAL cubed_sphere_rmax;

  // Swirl
  const int swirl_ncells_i;
  const int swirl_ncells_j;
  const int swirl_ncells_k;

  // Cake
  const CCTK_REAL
      cake_outer_boundary_radius; // The radius of the outer boundary

  const CCTK_REAL
      cake_inner_boundary_radius; /* Half the coordinate length of the
                                   * central cartesian cube's face
                                   */

  const int cake_cartesian_ncells_i; /* The number of cells in the x direction
                                      * of the central cartesian cube.
                                      */

  const int cake_cartesian_ncells_j; /* The number of cells in the y direction
                                      * of the central cartesian cube.
                                      */

  const int cake_cartesian_ncells_k; /* The number of cells in the z direction
                                      * of the central cartesian cube.
                                      */

  const int cake_angular_cells; /* The number of cells in the angular direction
                                 * of spherical patches.
                                 */

  const int cake_radial_cells; /* The number of cells in the radial direction of
                                * spherical patches.
                                */

  PatchTransformations();

  PatchTransformations(const PatchTransformations &) = default;
  PatchTransformations(PatchTransformations &&) = default;
  PatchTransformations &operator=(const PatchTransformations &) = default;
  PatchTransformations &operator=(PatchTransformations &&) = default;

  std::tuple<int, vec<CCTK_REAL, dim, UP> > (*global2local)(
      const PatchTransformations &pt, const vec<CCTK_REAL, dim, UP> &x) = 0;

  vec<CCTK_REAL, dim, UP> (*local2global)(const PatchTransformations &pt,
                                          int patch,
                                          const vec<CCTK_REAL, dim, UP> &a) = 0;

  // Calculating global derivatives d/dx from local derivatives d/da
  // requires the Jacobian da/dx, and also its derivative d^2/dx^2 for
  // second derivatives

  // da/dx[i,j] = da[i] / dx[j]
  std::tuple<vec<CCTK_REAL, dim, UP>, vec<vec<CCTK_REAL, dim, DN>, dim, UP> > (
      *dlocal_dglobal)(const PatchTransformations &pt, int patch,
                       const vec<CCTK_REAL, dim, UP> &a) = 0;

  // d^2a/dx^2[i,j,k] = d^2a[i] / dx^2[j,k]
  std::tuple<vec<CCTK_REAL, dim, UP>, vec<vec<CCTK_REAL, dim, DN>, dim, UP>,
             vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP> > (*d2local_dglobal2)(
      const PatchTransformations &pt, int patch,
      const vec<CCTK_REAL, dim, UP> &a) = 0;

  // Device functions mirroring the function above
  std::tuple<int, vec<CCTK_REAL, dim, UP> > (*global2local_device)(
      const PatchTransformations &pt, const vec<CCTK_REAL, dim, UP> &x) = 0;
  vec<CCTK_REAL, dim, UP> (*local2global_device)(
      const PatchTransformations &pt, int patch,
      const vec<CCTK_REAL, dim, UP> &a) = 0;
  std::tuple<vec<CCTK_REAL, dim, UP>, vec<vec<CCTK_REAL, dim, DN>, dim, UP> > (
      *dlocal_dglobal_device)(const PatchTransformations &pt, int patch,
                              const vec<CCTK_REAL, dim, UP> &a) = 0;
  std::tuple<vec<CCTK_REAL, dim, UP>, vec<vec<CCTK_REAL, dim, DN>, dim, UP>,
             vec<smat<CCTK_REAL, dim, DN, DN>, dim, UP> > (
      *d2local_dglobal2_device)(const PatchTransformations &pt, int patch,
                                const vec<CCTK_REAL, dim, UP> &a) = 0;
};

struct PatchSystem {
  std::vector<Patch> patches;
  int num_patches() const { return patches.size(); }

  PatchTransformations transformations;

  PatchSystem() {}
  PatchSystem(std::vector<Patch> patches, PatchTransformations transformations)
      : patches(std::move(patches)),
        transformations(std::move(transformations)) {}
};

////////////////////////////////////////////////////////////////////////////////

PatchSystem SetupCartesian();
PatchSystem SetupCubedSphere();
PatchSystem SetupSwirl();
PatchSystem SetupCake();

extern std::unique_ptr<PatchSystem> the_patch_system;

} // namespace MultiPatch

#endif // #ifndef MULTIPATCH_HXX
