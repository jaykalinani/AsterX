#include "weyl.hxx"

#include "weyl_vars.hxx"

namespace Weyl {

void gfs_t::calc_metric() const {
  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [layout0 = layout0, //
       gf_gamma0 = gf_gamma0, gf_alpha0 = gf_alpha0, gf_beta0 = gf_beta0,
       gf_K0 = gf_K0, gf_dtalpha0 = gf_dtalpha0, gf_dtbeta0 = gf_dtbeta0,
       gf_dgamma0 = gf_dgamma0, gf_dalpha0 = gf_dalpha0, gf_dbeta0 = gf_dbeta0,
       gf_dtK0 = gf_dtK0, gf_dt2alpha0 = gf_dt2alpha0,
       gf_dt2beta0 = gf_dt2beta0, gf_dK0 = gf_dK0, gf_ddtalpha0 = gf_ddtalpha0,
       gf_ddtbeta0 = gf_ddtbeta0, gf_ddgamma0 = gf_ddgamma0,
       gf_ddalpha0 = gf_ddalpha0, gf_ddbeta0 = gf_ddbeta0, //
       tile_g4 = tile_g4, tile_dg4 = tile_dg4,
       tile_ddg4 = tile_ddg4 //
  ] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index0(layout0, p.I);

        // Load and calculate

        const auto id3 = one<smat<int, 3> >()();

        const weyl_vars_metric<vreal> vars(
            gf_gamma0(mask, index0, id3), gf_alpha0(mask, index0, 1),
            gf_beta0(mask, index0), //
            gf_K0(mask, index0), gf_dtalpha0(mask, index0),
            gf_dtbeta0(mask, index0), //
            gf_dgamma0(mask, index0), gf_dalpha0(mask, index0),
            gf_dbeta0(mask, index0), //
            gf_dtK0(mask, index0), gf_dt2alpha0(mask, index0),
            gf_dt2beta0(mask, index0), //
            gf_dK0(mask, index0), gf_ddtalpha0(mask, index0),
            gf_ddtbeta0(mask, index0), //
            gf_ddgamma0(mask, index0), gf_ddalpha0(mask, index0),
            gf_ddbeta0(mask, index0));

        // Store

        tile_g4.store(mask, index0, vars.g);
        tile_dg4.store(mask, index0, vars.dg);
        tile_ddg4.store(mask, index0, vars.ddg);
      });
}

} // namespace Weyl
