#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cmath>

#include "prim2con.h"

namespace GRHydroToyGPU {
using namespace std;
using namespace Loop;

extern "C" void GRHydroToyGPU_Prim2Con_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_GRHydroToyGPU_Prim2Con_Initial;
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
  const GF3D2<const CCTK_REAL> gf_eps(gf_layout_cell, eps); 

  // Loop over the entire grid (0 to n-1 points in each direction)
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones, [=] CCTK_DEVICE CCTK_HOST(const PointDesc &p)
          CCTK_ATTRIBUTE_ALWAYS_INLINE {
  
          // Interpolate metric terms from vertices to center
	  metric g;
          g.gxx = 0.0;
          g.gxy = 0.0;
          g.gxz = 0.0;
          g.gyy = 0.0;
          g.gyz = 0.0;
          g.gzz = 0.0;

          for(int dk = 0 ; dk < 2 ; ++dk)
              for(int dj = 0 ; dj < 2 ; ++dj)
                  for(int di = 0 ; di < 2 ; ++di) {
                      g.gxx += gf_gxx(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      g.gxy += gf_gxy(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      g.gxz += gf_gxz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      g.gyy += gf_gyy(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      g.gyz += gf_gyz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                      g.gzz += gf_gzz(p.I + p.DI[0]*di + p.DI[1]*dj + p.DI[2]*dk);
                  }

          g.gxx *= 0.125;
          g.gxy *= 0.125;
          g.gxz *= 0.125;
          g.gyy *= 0.125;
          g.gyz *= 0.125;
          g.gzz *= 0.125;
	  
 
	  prim pv;
	  pv.rho   = gf_rho(p.I);
	  pv.velx  = gf_velx(p.I);
	  pv.vely  = gf_vely(p.I);
	  pv.velz  = gf_velz(p.I);
	  pv.eps   = gf_eps(p.I);
	  pv.press = gf_press(p.I);
          
	  cons cv;
	  prim2con(g,pv,cv);

	  gf_dens(p.I) = cv.dens;
	  gf_momx(p.I) = cv.momx;
	  gf_momy(p.I) = cv.momy;
	  gf_momz(p.I) = cv.momz;
	  gf_tau(p.I)  = cv.tau;

        });
    }
} // namespace GRHydroToyGPU
