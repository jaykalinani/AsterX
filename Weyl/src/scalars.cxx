#include "weyl.hxx"

#include "weyl_vars.hxx"

namespace Weyl {

void gfs_t::calc_scalars() const {
  typedef simd<CCTK_REAL> vreal;
  typedef simdl<CCTK_REAL> vbool;
  constexpr std::size_t vsize = std::tuple_size_v<vreal>;

  const CCTK_REAL time = cctkGH->cctk_time;

  grid.loop_int_device<0, 0, 0, vsize>(
      grid.nghostzones, [layout0 = layout0, layout5 = layout5, //
                         time,                                 //
                         tile_g4 = tile_g4, tile_R4 = tile_R4,
                         tile_C4 = tile_C4, //
                         gf_Psi0re5 = gf_Psi0re5, gf_Psi0im5 = gf_Psi0im5,
                         gf_Psi1re5 = gf_Psi1re5, gf_Psi1im5 = gf_Psi1im5,
                         gf_Psi2re5 = gf_Psi2re5, gf_Psi2im5 = gf_Psi2im5,
                         gf_Psi3re5 = gf_Psi3re5, gf_Psi3im5 = gf_Psi3im5,
                         gf_Psi4re5 = gf_Psi4re5,
                         gf_Psi4im5 = gf_Psi4im5 //
  ] ARITH_DEVICE(const PointDesc &p) ARITH_INLINE {
        const vbool mask = mask_for_loop_tail<vbool>(p.i, p.imax);
        const GF3D5index index0(layout0, p.I);
        const GF3D5index index5(layout5, p.I);

        // Load and calculate

        constexpr smat<int, 4> id4{-1, 0, 0, 0, //
                                   +1, 0, 0,    //
                                   +1, 0,       //
                                   +1};

        const vec<vreal, 3> coord3(
            [&](int d) { return p.X[d] + iota<vreal>() * p.DX[d]; });
        const vec<vreal, 4> coord4{time, coord3(0), coord3(1), coord3(2)};

        const weyl_vars_scalars<vreal> vars(coord4, //
                                            tile_g4(mask, index0, id4),
                                            tile_R4(mask, index0),
                                            tile_C4(mask, index0));

        // Store

        gf_Psi0re5.store(mask, index5, real(vars.Psi0));
        gf_Psi0im5.store(mask, index5, imag(vars.Psi0));
        gf_Psi1re5.store(mask, index5, real(vars.Psi1));
        gf_Psi1im5.store(mask, index5, imag(vars.Psi1));
        gf_Psi2re5.store(mask, index5, real(vars.Psi2));
        gf_Psi2im5.store(mask, index5, imag(vars.Psi2));
        gf_Psi3re5.store(mask, index5, real(vars.Psi3));
        gf_Psi3im5.store(mask, index5, imag(vars.Psi3));
        gf_Psi4re5.store(mask, index5, real(vars.Psi4));
        gf_Psi4im5.store(mask, index5, imag(vars.Psi4));
      });
}

} // namespace Weyl
