#ifndef WEYL_HXX
#define WEYL_HXX

#include <loop_device.hxx>
#include <mat.hxx>
#include <rten.hxx>
#include <vec.hxx>
#include <vect.hxx>

#include <cctk.h>

#include <array>
#include <type_traits>

namespace Weyl {
using namespace Arith;
using namespace Loop;

struct gfs_t {

  const cGH *restrict cctkGH;

  std::array<int, dim> indextype;
  std::array<int, dim> nghostzones;
  GridDescBaseDevice grid;
  vect<int, dim> imin, imax;

  // Suffix 1: with ghost zones,
  // suffix 5: with ghost zones,
  // suffix 0: without ghost zones
  GF3D2layout layout1;
  GF3D5layout layout5;
  GF3D5layout layout0;

  // Input grid functions

  GF3D2<const CCTK_REAL> gf_alpha1;
  vec<GF3D2<const CCTK_REAL>, 3> gf_beta1;
  smat<GF3D2<const CCTK_REAL>, 3> gf_gamma1;
  smat<GF3D2<const CCTK_REAL>, 3> gf_K1;

  GF3D2<const CCTK_REAL> gf_dtalpha1;
  vec<GF3D2<const CCTK_REAL>, 3> gf_dtbeta1;
  smat<GF3D2<const CCTK_REAL>, 3> gf_dtK1;

  GF3D2<const CCTK_REAL> gf_dt2alpha1;
  vec<GF3D2<const CCTK_REAL>, 3> gf_dt2beta1;

  // Output grid functions

  GF3D5<CCTK_REAL> gf_Psi0re5;
  GF3D5<CCTK_REAL> gf_Psi0im5;
  GF3D5<CCTK_REAL> gf_Psi1re5;
  GF3D5<CCTK_REAL> gf_Psi1im5;
  GF3D5<CCTK_REAL> gf_Psi2re5;
  GF3D5<CCTK_REAL> gf_Psi2im5;
  GF3D5<CCTK_REAL> gf_Psi3re5;
  GF3D5<CCTK_REAL> gf_Psi3im5;
  GF3D5<CCTK_REAL> gf_Psi4re5;
  GF3D5<CCTK_REAL> gf_Psi4im5;

private:
  // Temporary variables

  int nvars;
  mutable int ivar;
  GF3D5vector<CCTK_REAL> vars;

  template <typename F, typename R = std::result_of_t<F()> >
  static auto make_vec(const F &f) {
    return vec<R, 3>([&](int) { return f(); });
  }
  template <typename F, typename R = std::result_of_t<F()> >
  static auto make_mat(const F &f) {
    return smat<R, 3>([&](int, int) { return f(); });
  }
  template <typename F, typename R = std::result_of_t<F()> >
  static auto make_vec4(const F &f) {
    return vec<R, 4>([&](int) { return f(); });
  }
  template <typename F, typename R = std::result_of_t<F()> >
  static auto make_mat4(const F &f) {
    return smat<R, 4>([&](int, int) { return f(); });
  }
  template <typename F, typename R = std::result_of_t<F()> >
  static auto make_rten4(const F &f) {
    return rten<R, 4>([&](int, int, int, int) { return f(); });
  }

  auto make_gf() const { return GF3D5<CCTK_REAL>(vars(ivar++)); }
  auto make_vec_gf() const {
    return make_vec([&]() { return make_gf(); });
  }
  auto make_mat_gf() const {
    return make_mat([&]() { return make_gf(); });
  }
  auto make_vec_vec_gf() const {
    return make_vec([&]() { return make_vec_gf(); });
  }
  auto make_vec_mat_gf() const {
    return make_vec([&]() { return make_mat_gf(); });
  }
  auto make_mat_vec_gf() const {
    return make_mat([&]() { return make_vec_gf(); });
  }
  auto make_mat_mat_gf() const {
    return make_mat([&]() { return make_mat_gf(); });
  }

  auto make_vec4_gf() const {
    return make_vec4([&]() { return make_gf(); });
  }
  auto make_mat4_gf() const {
    return make_mat4([&]() { return make_gf(); });
  }
  auto make_rten4_gf() const {
    return make_rten4([&]() { return make_gf(); });
  }
  auto make_vec4_mat4_gf() const {
    return make_vec4([&]() { return make_mat4_gf(); });
  }
  auto make_mat4_vec4_gf() const {
    return make_mat4([&]() { return make_vec4_gf(); });
  }
  auto make_mat4_mat4_gf() const {
    return make_mat4([&]() { return make_mat4_gf(); });
  }

public:
  // Input variables: ADM variables

  GF3D5<CCTK_REAL> gf_alpha0;
  vec<GF3D5<CCTK_REAL>, 3> gf_dalpha0;
  smat<GF3D5<CCTK_REAL>, 3> gf_ddalpha0;

  vec<GF3D5<CCTK_REAL>, 3> gf_beta0;
  vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dbeta0;
  vec<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddbeta0;

  smat<GF3D5<CCTK_REAL>, 3> gf_gamma0;
  smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dgamma0;
  smat<smat<GF3D5<CCTK_REAL>, 3>, 3> gf_ddgamma0;

  smat<GF3D5<CCTK_REAL>, 3> gf_K0;
  smat<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_dK0;

  GF3D5<CCTK_REAL> gf_dtalpha0;
  vec<GF3D5<CCTK_REAL>, 3> gf_ddtalpha0;

  vec<GF3D5<CCTK_REAL>, 3> gf_dtbeta0;
  vec<vec<GF3D5<CCTK_REAL>, 3>, 3> gf_ddtbeta0;

  smat<GF3D5<CCTK_REAL>, 3> gf_dtK0;

  GF3D5<CCTK_REAL> gf_dt2alpha0;

  vec<GF3D5<CCTK_REAL>, 3> gf_dt2beta0;

  // Intermediate variables: 4-metric

  smat<GF3D5<CCTK_REAL>, 4> tile_g4;
  smat<vec<GF3D5<CCTK_REAL>, 4>, 4> tile_dg4;
  smat<smat<GF3D5<CCTK_REAL>, 4>, 4> tile_ddg4;

  // Intermediate variables: 4-curvature

  vec<smat<GF3D5<CCTK_REAL>, 4>, 4> tile_Gamma4;
  smat<GF3D5<CCTK_REAL>, 4> tile_R4;
  rten<GF3D5<CCTK_REAL>, 4> tile_C4;

  //

  gfs_t() = delete;
  gfs_t(const gfs_t &) = delete;
  gfs_t(gfs_t &&) = delete;
  gfs_t &operator=(const gfs_t &) = delete;
  gfs_t &operator=(gfs_t &&) = delete;

  gfs_t(const cGH *cctkGH);

  void calc_metric() const;
  void calc_curvature() const;
  void calc_scalars() const;
};

} // namespace Weyl

#endif // #ifndef WEYL_HXX
