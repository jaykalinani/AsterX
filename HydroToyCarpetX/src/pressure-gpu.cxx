//#include <loop.hxx>
//#include <vectors.h>
//
#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>
#include <cctk_Arguments.h>
#include <loop_device.hxx>

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>

namespace HydroToyCarpetX {
using namespace std;

constexpr int dim = 3;

namespace {
template <typename T> T pow2(T x) { return x * x; }
} // namespace

extern "C" void HydroToyCarpetX_Pressure(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_HydroToyCarpetX_Pressure;
  DECLARE_CCTK_PARAMETERS;

  // regular, cell centred variables

  // assert((imax[0] - imin[0]) % VS == 0);

  // Equation of state: p = (gamma - 1) e

  // vel^j = delta^j_i mom_i / rho

//  printf("Before loop\n");

  // Determine loop extent
  const array<int, dim> nghostzones{cctkGH->cctk_nghostzones[0],
                                    cctkGH->cctk_nghostzones[1],
                                    cctkGH->cctk_nghostzones[2]};
  const Loop::GridDescBaseDevice griddesc(cctkGH);
  griddesc.loop_all_device<1, 1, 1>(
      nghostzones, [=] CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_HOST CCTK_DEVICE(
                       const Loop::PointDesc &p) {

        CCTK_REAL rho_inv = 1.0/rho[p.idx];
//        printf("rho = %g\n", rho[p.idx]);
//        printf("etot = %g\n", etot[p.idx]);

        CCTK_REAL ekin = 0.5 * rho_inv * sqrt(momx[p.idx]*momx[p.idx] + momy[p.idx]*momy[p.idx] + momz[p.idx]*momz[p.idx]);
        eint[p.idx] = etot[p.idx] - ekin;

        press[p.idx] = (gamma - 1) * eint[p.idx];

//        printf("p = %g\n", press[p.idx]);
        velx[p.idx] = rho_inv * momx[p.idx];
        vely[p.idx] = rho_inv * momy[p.idx];
        velz[p.idx] = rho_inv * momz[p.idx];
     }); 
  }
} // namespace HydroToyCarpetX
