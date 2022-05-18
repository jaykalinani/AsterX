#include <loop_device.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>
#include <random>

//initial setup adapted from Springel+2010

namespace KHInitial {
using namespace std;
using namespace Loop;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T pow2(T x) {
  return x * x;
}

extern "C" void KHInitial_Initialize(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_KHInitial_Initialize;
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

  if (CCTK_EQUALS(initial_hydro, "KHI")) {

    grid.loop_all_device<1, 1, 1>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

	  //set velocities in opposite directions in x along the slip surface 
	  //set slightly different densities along the slip surface
	  //location of slip surfaces hardcoded at (abs(p.y) = 0.25)
          if (abs(p.y) >= 0.25) {
            gf_rho(p.I) = rhoUp;
            gf_velx(p.I) = vxUp;
          } else {
            gf_rho(p.I) = rhoLow;
            gf_velx(p.I) = vxLow;
          }
          //excite the instability by peturbing v^y
	  gf_vely(p.I) = w0*sin(4*M_PI*p.x) * ( exp(-pow2(p.y-0.25)/(2 * pow2(sigma))) + exp(-pow2(p.y+0.25)/(2*pow2(sigma))) );
	  gf_velz(p.I) = 0.0;
          
	  //set constant initial pressure throughout the domain
	  gf_press(p.I) = p_val;

	  // TODO: compute eps using EOS driver
          // for now, using ideal gas EOS
          gf_eps(p.I) = gf_press(p.I) / (gf_rho(p.I) * (gamma - 1));
        });

  } else {
    CCTK_ERROR("Incorrect test case!");
  }
}

} // namespace KHInitial
