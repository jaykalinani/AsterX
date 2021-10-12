#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

////////////////////////////////////////////////////////////////////////////////

extern "C" void GRHydroToyGPU_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Initialize;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  const GF3D2layout gf_layout_cell(cctkGH, cell_centred);

  const GF3D2<CCTK_REAL> gf_rho(gf_layout_cell, rho);
  const GF3D2<CCTK_REAL> gf_velx(gf_layout_cell, velx);
  const GF3D2<CCTK_REAL> gf_vely(gf_layout_cell, vely);
  const GF3D2<CCTK_REAL> gf_velz(gf_layout_cell, velz);
  const GF3D2<CCTK_REAL> gf_press(gf_layout_cell, press);
  const GF3D2<CCTK_REAL> gf_eps(gf_layout_cell, eps);

  if (CCTK_EQUALS(setup, "equilibrium")) {

    grid.loop_int_device<1, 1, 1>(
	grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
            CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gf_rho(p.I) = 1.0;
          gf_velx(p.I) = 0.0;
          gf_vely(p.I) = 0.0;
          gf_velz(p.I) = 0.0;
          gf_press(p.I) = 1.0;
          //TODO: compute eps using EOS driver
	  //for now, using ideal gas EOS
	  gf_eps(p.I) = gf_press(p.I)/(gf_rho(p.I)*(gamma-1) );
    });

  } else if (CCTK_EQUALS(setup, "sound wave")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                              const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          gf_rho(p.I) = 1.0;
          gf_velx(p.I) = 0.0 + amplitude * sin(M_PI * p.x);
          gf_vely(p.I) = 0.0;
          gf_velz(p.I) = 0.0;
          gf_press(p.I) = 1.0; // should add kinetic energy here
	  //TODO: compute eps using EOS driver
          //for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I)/(gf_rho(p.I)*(gamma-1) );
        });

  } else if (CCTK_EQUALS(setup, "shock tube")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
	    const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
                if (p.x <= 0.0) {
                    gf_rho(p.I) = 2.0;
                    gf_velx(p.I) = 0.0;
                    gf_vely(p.I) = 0.0;
                    gf_velz(p.I) = 0.0;
                    gf_press(p.I) = 2.0;
                } else {
                    gf_rho(p.I) = 1.0;
                    gf_velx(p.I) = 0.0;
                    gf_vely(p.I) = 0.0;
                    gf_velz(p.I) = 0.0;
                    gf_press(p.I) = 1.0;
                } 

          //TODO: compute eps using EOS driver
          //for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I)/(gf_rho(p.I)*(gamma-1) );
		
      });

  } else if (CCTK_EQUALS(setup, "balsara1")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
            const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
                if (p.x <= 0.0) {
                    gf_rho(p.I) = 1.0;
                    gf_velx(p.I) = 0.0;
                    gf_vely(p.I) = 0.0;
                    gf_velz(p.I) = 0.0;
                    gf_press(p.I) = 1.0;
                } else {
                    gf_rho(p.I) = 0.125;
                    gf_velx(p.I) = 0.0;
                    gf_vely(p.I) = 0.0;
                    gf_velz(p.I) = 0.0;
                    gf_press(p.I) = 0.1;
                }

          //TODO: compute eps using EOS driver
          //for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I)/(gf_rho(p.I)*(gamma-1) );

      });

  } else if (CCTK_EQUALS(setup, "spherical shock")) {

    grid.loop_int_device<1, 1, 1>(
        grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(
                              const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
          CCTK_REAL r2 = pow(p.x, 2) + pow(p.y, 2) + pow(p.z, 2);
          if (r2 <= pow(shock_radius, 2)) {
            gf_rho(p.I) = 2.0;
            gf_velx(p.I) = 0.0;
            gf_vely(p.I) = 0.0;
            gf_velz(p.I) = 0.0;
            gf_press(p.I) = 2.0;
          } else {
            gf_rho(p.I) = 1.0;
            gf_velx(p.I) = 0.0;
            gf_vely(p.I) = 0.0;
            gf_velz(p.I) = 0.0;
            gf_press(p.I) = 1.0;
          }

	  //TODO: compute eps using EOS driver
          //for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I)/(gf_rho(p.I)*(gamma-1) );

        });

  } else {
    CCTK_ERROR("Internal error");
  }
}

} // namespace GRHydroToyGPU
