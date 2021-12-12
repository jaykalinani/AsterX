#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <reprimand/con2prim_imhd.h>
#include <reprimand/eos_idealgas.h>

#include <cmath>
#include <iostream> //TODO

namespace HydroToyGPU {
using namespace std;
using namespace EOS_Toolkit;
using namespace Loop;

namespace {
template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
}
} // namespace

extern "C" void HydroToyGPU_Pressure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyGPU_Pressure;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout(cctkGH, cell_centred);

  const GF3D2<const CCTK_REAL> gf_rho(gf_layout, rho);
  const GF3D2<const CCTK_REAL> gf_momx(gf_layout, momx);
  const GF3D2<const CCTK_REAL> gf_momy(gf_layout, momy);
  const GF3D2<const CCTK_REAL> gf_momz(gf_layout, momz);
  const GF3D2<const CCTK_REAL> gf_etot(gf_layout, etot);

  const GF3D2<CCTK_REAL> gf_press(gf_layout, press);
  const GF3D2<CCTK_REAL> gf_velx(gf_layout, velx);
  const GF3D2<CCTK_REAL> gf_vely(gf_layout, vely);
  const GF3D2<CCTK_REAL> gf_velz(gf_layout, velz);
  const GF3D2<CCTK_REAL> gf_eint(gf_layout, eint);

  if (CCTK_EQUALS(con2prim_method, "direct")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL rho_inv = 1.0 / (gf_rho(p.I) + 1.0e-20);

          CCTK_REAL ekin =
              0.5 * rho_inv *
              (pow2(gf_momx(p.I)) + pow2(gf_momy(p.I)) + pow2(gf_momz(p.I)));
          gf_eint(p.I) = gf_etot(p.I) - ekin;

          gf_press(p.I) = (gamma - 1) * gf_eint(p.I);

          gf_velx(p.I) = rho_inv * gf_momx(p.I);
          gf_vely(p.I) = rho_inv * gf_momy(p.I);
          gf_velz(p.I) = rho_inv * gf_momz(p.I);
        });

  } else if (CCTK_EQUALS(con2prim_method, "reprimand")) {

    // Define the EOS
    const CCTK_REAL max_eps = 11.0;
    const CCTK_REAL max_rho = 1.0e6;
    const CCTK_REAL adiab_ind = 1.0; // TODO: use gamma
    const auto eos = make_eos_idealgas(adiab_ind, max_eps, max_rho);

    // Set up atmosphere
    const CCTK_REAL atmo_rho = 1.0e-20;
    const CCTK_REAL atmo_eps = 0.1;
    const CCTK_REAL atmo_ye = 0.5;
    const CCTK_REAL atmo_cut = atmo_rho * 1.01;
    const CCTK_REAL atmo_p =
        eos.at_rho_eps_ye(atmo_rho, atmo_eps, atmo_ye).press();
    const atmosphere atmo{atmo_rho, atmo_eps, atmo_ye, atmo_p, atmo_cut};

    // Primitive recovery
    const CCTK_REAL rho_strict = 1.0e-11;
    const bool ye_lenient = false;
    const int max_iter = 30;
    const CCTK_REAL c2p_acc = 1.0e-8;
    const CCTK_REAL max_b = 10.0;
    const CCTK_REAL max_z = 1.0e3;
    const con2prim_mhd cv2pv(eos, rho_strict, ye_lenient, max_z, max_b, atmo,
                             c2p_acc, max_iter);

    grid.loop_all<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_HOST(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          sm_metric3 g;
          g.minkowski();

          const CCTK_REAL ekin =
              0.5 / gf_rho(p.I) *
              (pow2(gf_momx(p.I)) + pow2(gf_momy(p.I)) + pow2(gf_momz(p.I)));
          const CCTK_REAL eint = gf_etot(p.I) - ekin;
          cons_vars_mhd conserved{gf_rho(p.I),
                                  eint / gf_rho(p.I),
                                  0.5 * gf_rho(p.I),
                                  {gf_momx(p.I), gf_momy(p.I), gf_momz(p.I)},
                                  {0.0, 0.0, 0.0}};
          // std::cout << "dens=" << conserved.dens << "\n"
          //           << "tau=" << conserved.tau << "\n"
          //           << "tracer_ye=" << conserved.tracer_ye << "\n"
          //           << "scon=[" << conserved.scon(0) << "," <<
          //           conserved.scon(1)
          //           << "," << conserved.scon(2) << "]\n"
          //           << "bcons=[" << conserved.bcons(0) << ","
          //           << conserved.bcons(1) << "," << conserved.bcons(2) <<
          //           "]\n";

          prim_vars_mhd primitives;
          con2prim_mhd::report rep;
          cv2pv(primitives, conserved, g, rep);

          assert(!rep.failed());
          if (rep.failed())
            CCTK_ERROR(rep.debug_message().c_str());

          // gf_rho(p.I) = primitives.rho;
          gf_velx(p.I) = primitives.vel(0);
          gf_vely(p.I) = primitives.vel(1);
          gf_velz(p.I) = primitives.vel(2);
          gf_press(p.I) = primitives.press;
          gf_eint(p.I) = primitives.rho * primitives.eps;

          if (rep.adjust_cons) {
            assert(0);
            //   gf_rho(p.I) = conserved.dens;
            //   gf_momx(p.I) = conserved.scon(0);
            //   gf_momy(p.I) = conserved.scon(1);
            //   gf_momz(p.I) = conserved.scon(2);
            //   gf_etot(p.I) = conserved.tau;
          }
        });

  } else {
    CCTK_ERROR("internal error");
  }
}

} // namespace HydroToyGPU
