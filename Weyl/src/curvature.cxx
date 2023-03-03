#include "weyl.hxx"

#include "weyl_vars.hxx"

namespace Weyl {

void gfs_t::calc_curvature() const {
  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones,
      [layout0 = layout0,                                             //
       tile_g4 = tile_g4, tile_dg4 = tile_dg4, tile_ddg4 = tile_ddg4, //
       // tile_Gamma4 = tile_Gamma4, tile_R4 = tile_R4,
       tile_C4 = tile_C4] //
      ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index0(layout0, p.I);

        // Load and calculate

        constexpr smat<int, 4> id4{-1, 0, 0, 0, //
                                   +1, 0, 0,    //
                                   +1, 0,       //
                                   +1};

        const weyl_vars_curvature<vreal> vars(tile_g4(mask, index0, id4),
                                              tile_dg4(mask, index0),
                                              tile_ddg4(mask, index0));

        // Store

        // tile_Gamma4.store(mask, index0, vars.Gamma);
        // tile_R4.store(mask, index0, vars.R);

        // tile_C4.store(mask, index0, vars.C);
        tile_C4.get(0, 1, 0, 1).store(mask, index0, vars.C(0, 1, 0, 1));
        tile_C4.get(0, 1, 0, 2).store(mask, index0, vars.C(0, 1, 0, 2));
        tile_C4.get(0, 1, 0, 3).store(mask, index0, vars.C(0, 1, 0, 3));
        tile_C4.get(0, 1, 1, 2).store(mask, index0, vars.C(0, 1, 1, 2));
        tile_C4.get(0, 1, 1, 3).store(mask, index0, vars.C(0, 1, 1, 3));
        tile_C4.get(0, 1, 2, 3).store(mask, index0, vars.C(0, 1, 2, 3));

        tile_C4.get(0, 2, 0, 2).store(mask, index0, vars.C(0, 2, 0, 2));
        tile_C4.get(0, 2, 0, 3).store(mask, index0, vars.C(0, 2, 0, 3));
        tile_C4.get(0, 2, 1, 2).store(mask, index0, vars.C(0, 2, 1, 2));
        tile_C4.get(0, 2, 1, 3).store(mask, index0, vars.C(0, 2, 1, 3));
        tile_C4.get(0, 2, 2, 3).store(mask, index0, vars.C(0, 2, 2, 3));

        tile_C4.get(0, 3, 0, 3).store(mask, index0, vars.C(0, 3, 0, 3));
        tile_C4.get(0, 3, 1, 2).store(mask, index0, vars.C(0, 3, 1, 2));
        tile_C4.get(0, 3, 1, 3).store(mask, index0, vars.C(0, 3, 1, 3));
        tile_C4.get(0, 3, 2, 3).store(mask, index0, vars.C(0, 3, 2, 3));

        tile_C4.get(1, 2, 1, 2).store(mask, index0, vars.C(1, 2, 1, 2));
        tile_C4.get(1, 2, 1, 3).store(mask, index0, vars.C(1, 2, 1, 3));
        tile_C4.get(1, 2, 2, 3).store(mask, index0, vars.C(1, 2, 2, 3));

        tile_C4.get(1, 3, 1, 3).store(mask, index0, vars.C(1, 3, 1, 3));
        tile_C4.get(1, 3, 2, 3).store(mask, index0, vars.C(1, 3, 2, 3));

        tile_C4.get(2, 3, 2, 3).store(mask, index0, vars.C(2, 3, 2, 3));
      });
}

} // namespace Weyl
