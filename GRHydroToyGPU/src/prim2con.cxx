#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

extern "C" void GRHydroToyGPU_Prim2Con(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Prim2Con;
  DECLARE_CCTK_PARAMETERS;

  const GridDescBaseDevice grid(cctkGH);
  constexpr array<int, dim> cell_centred = {1, 1, 1};
  constexpr array<int, dim> vertex_centred = {0, 0, 0};
  const GF3D2layout gf_layout_cell(cctkGH, cell_centred);
  const GF3D2layout gf_layout_vertex(cctkGH, vertex_centred);

  const GF3D2<const CCTK_REAL> gf_gxx(gf_layout_vertex, gxx);
  const GF3D2<const CCTK_REAL> gf_gxy(gf_layout_vertex, gxy);
  const GF3D2<const CCTK_REAL> gf_gxz(gf_layout_vertex, gxz);
  const GF3D2<const CCTK_REAL> gf_gyy(gf_layout_vertex, gyy);
  const GF3D2<const CCTK_REAL> gf_gyz(gf_layout_vertex, gyz);
  const GF3D2<const CCTK_REAL> gf_gzz(gf_layout_vertex, gzz);

  GF3D2<CCTK_REAL> gf_dens(gf_layout_cell, dens);
  GF3D2<CCTK_REAL> gf_momx(gf_layout_cell, momx);
  GF3D2<CCTK_REAL> gf_momy(gf_layout_cell, momy);
  GF3D2<CCTK_REAL> gf_momz(gf_layout_cell, momz);
  GF3D2<CCTK_REAL> gf_tau(gf_layout_cell, tau);

  const GF3D2<const CCTK_REAL> gf_rho(gf_layout_cell, rho);
  const GF3D2<const CCTK_REAL> gf_velx(gf_layout_cell, velx);
  const GF3D2<const CCTK_REAL> gf_vely(gf_layout_cell, vely);
  const GF3D2<const CCTK_REAL> gf_velz(gf_layout_cell, velz);
  const GF3D2<const CCTK_REAL> gf_press(gf_layout_cell, press);
  GF3D2<CCTK_REAL> gf_eps(gf_layout_cell, eps); 

  // Loop over the entire grid (0 to n-1 points in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
  
          // Interpolate metric terms from vertices to center
          CCTK_REAL gxx_avg = 0.0;
          CCTK_REAL gxy_avg = 0.0;
          CCTK_REAL gxz_avg = 0.0;
          CCTK_REAL gyy_avg = 0.0;
          CCTK_REAL gyz_avg = 0.0;
          CCTK_REAL gzz_avg = 0.0;

          for(int dk = 0 ; dk < 2 ; ++dk)
              for(int dj = 0 ; dj < 2 ; ++dj)
                  for(int di = 0 ; di < 2 ; ++di) {
                      gxx_avg += gf_gxx(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      gxy_avg += gf_gxy(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      gxz_avg += gf_gxz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      gyy_avg += gf_gyy(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      gyz_avg += gf_gyz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      gzz_avg += gf_gzz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                  }

          gxx_avg *= 0.125;
          gxy_avg *= 0.125;
          gxz_avg *= 0.125;
          gyy_avg *= 0.125;
          gyz_avg *= 0.125;
          gzz_avg *= 0.125;
  
          //determinant of spatial metric
          const CCTK_REAL detg = -gxz_avg*gxz_avg*gyy_avg + 2.0*gxy_avg*gxz_avg*gyz_avg - gxx_avg*gyz_avg*gyz_avg
                                 - gxy_avg*gxy_avg*gzz_avg + gxx_avg*gyy_avg*gzz_avg;
          const CCTK_REAL sqrt_detg = sqrt(detg);
          
	  //TODO: compute specific internal energy based on user-specified EOS  
          //currently, computing eps for classical ideal gas
	
          gf_eps(p.I) = gf_press(p.I)/(gf_rho(p.I)*(gamma-1));

	  //v_j
          const CCTK_REAL vlowx = gxx_avg*gf_velx(p.I) + gxy_avg*gf_vely(p.I) + gxz_avg*gf_velz(p.I);
	  const CCTK_REAL vlowy = gxy_avg*gf_velx(p.I) + gyy_avg*gf_vely(p.I) + gyz_avg*gf_velz(p.I);
	  const CCTK_REAL vlowz = gxz_avg*gf_velx(p.I) + gyz_avg*gf_vely(p.I) + gzz_avg*gf_velz(p.I);

	  //w_lorentz
	  const CCTK_REAL w_lorentz = 1.0 / sqrt(1 - (gxx_avg*gf_velx(p.I)*gf_velx(p.I) 
				     + gyy_avg*gf_vely(p.I)*gf_vely(p.I)
			             + gzz_avg*gf_velz(p.I)*gf_velz(p.I) + 2*gxy_avg*gf_velx(p.I)*gf_vely(p.I)
			             + 2*gxz_avg*gf_velx(p.I)*gf_velz(p.I) + 2*gyz_avg*gf_vely(p.I)*gf_velz(p.I) ));


	  //computing conservatives from primitives
          gf_dens(p.I) = sqrt_detg * gf_rho(p.I)* w_lorentz;

          gf_momx(p.I) = sqrt_detg * gf_rho(p.I) * w_lorentz * (1 + gf_eps(p.I) 
				+ gf_press(p.I)/gf_rho(p.I)) * vlowx;
      
	  gf_momy(p.I) = sqrt_detg * gf_rho(p.I) * w_lorentz * (1 + gf_eps(p.I)    
                                + gf_press(p.I)/gf_rho(p.I)) * vlowy;

	  gf_momz(p.I) = sqrt_detg * gf_rho(p.I) * w_lorentz * (1 + gf_eps(p.I)
                                + gf_press(p.I)/gf_rho(p.I)) * vlowz;
  
	  gf_tau(p.I) = sqrt_detg * gf_rho(p.I) * w_lorentz * ( (1 + gf_eps(p.I)
                                + gf_press(p.I)/gf_rho(p.I)) * w_lorentz - 1) - gf_press(p.I);

        });
    }
} // namespace GRHydroToyGPU
