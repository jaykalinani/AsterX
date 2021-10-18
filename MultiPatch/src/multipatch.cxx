#include "multipatch.hxx"

#include <loop.hxx>

#include <cctk_Parameters.h>

#include <array>
#include <cmath>

namespace MultiPatch {

std::unique_ptr<PatchSystem> the_patch_system;

template <typename T, int D, dnup_t dnup1, dnup_t dnup2>
static T det(const vec<vec<T, D, dnup1>, D, dnup2> &A) {
  return A(0)(0) * (A(1)(1) * A(2)(2) - A(1)(2) * A(2)(1)) +
         A(0)(1) * (A(1)(2) * A(2)(0) - A(1)(0) * A(2)(2)) +
         A(0)(2) * (A(1)(0) * A(2)(1) - A(1)(1) * A(2)(0));
}

extern "C" int MultiPatch_Setup() {
  DECLARE_CCTK_PARAMETERS;

  if (CCTK_EQUALS(patch_system, "Cartesian")) {
    the_patch_system = SetupCartesian();
  } else if (CCTK_EQUALS(patch_system, "Cubed sphere")) {
    the_patch_system = SetupCubedSphere();
  } else if (CCTK_EQUALS(patch_system, "Swirl")) {
    the_patch_system = SetupSwirl();
  } else {
    CCTK_VERROR("Unknown patch system \"%s\"", patch_system);
  }

  return 0;
}

extern "C" void MultiPatch_Coordinates_Setup(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_MultiPatch_Coordinates_Setup;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GridDescBase grid(cctkGH);

  const std::array<int, dim> indextype_vc = {0, 0, 0};
  const std::array<int, dim> indextype_cc = {1, 1, 1};
  const Loop::GF3D2layout layout_vc(cctkGH, indextype_vc);
  const Loop::GF3D2layout layout_cc(cctkGH, indextype_cc);

  const Loop::GF3D2<CCTK_REAL> gf_vcoordx(layout_vc, vcoordx);
  const Loop::GF3D2<CCTK_REAL> gf_vcoordy(layout_vc, vcoordy);
  const Loop::GF3D2<CCTK_REAL> gf_vcoordz(layout_vc, vcoordz);

  const Loop::GF3D2<CCTK_REAL> gf_ccoordx(layout_cc, ccoordx);
  const Loop::GF3D2<CCTK_REAL> gf_ccoordy(layout_cc, ccoordy);
  const Loop::GF3D2<CCTK_REAL> gf_ccoordz(layout_cc, ccoordz);

  const Loop::GF3D2<CCTK_REAL> gf_cvol(layout_cc, cvol);

  grid.loop_all<0, 0, 0>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const Loop::PointDesc &p) ARITH_INLINE {
        const Loop::GF3D2index index(layout_vc, p.I);
        const vec<CCTK_REAL, dim, UP> a = {p.x, p.y, p.z};
        const vec<CCTK_REAL, dim, UP> x =
            the_patch_system->local2global(cctk_patch, a);
        gf_vcoordx(index) = x(0);
        gf_vcoordy(index) = x(1);
        gf_vcoordz(index) = x(2);
      });

  grid.loop_all<1, 1, 1>(
      grid.nghostzones,
      [=] ARITH_DEVICE ARITH_HOST(const Loop::PointDesc &p) ARITH_INLINE {
        const Loop::GF3D2index index(layout_cc, p.I);
        const vec<CCTK_REAL, dim, UP> a = {p.x, p.y, p.z};
        const std::tuple<vec<CCTK_REAL, dim, UP>,
                         vec<vec<CCTK_REAL, dim, DN>, dim, UP> >
            x_dadx = the_patch_system->dlocal_dglobal(cctk_patch, a);
        const vec<CCTK_REAL, dim, UP> &x = std::get<0>(x_dadx);
        const vec<vec<CCTK_REAL, dim, DN>, dim, UP> &dadx = std::get<1>(x_dadx);
        const CCTK_REAL det_dadx = det(dadx);
        using std::sqrt;
        const CCTK_REAL vol = (p.dx * p.dy * p.dz) * sqrt(det_dadx);
        gf_ccoordx(index) = x(0);
        gf_ccoordy(index) = x(1);
        gf_ccoordz(index) = x(2);
        gf_cvol(index) = vol;
      });
}

} // namespace MultiPatch
