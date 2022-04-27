#include "prim2con.hxx"

#include <fixmath.hxx>
#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

namespace AsterX {
using namespace Loop;

extern "C" void AsterX_Prim2Con_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterX_Prim2Con_Initial;
  DECLARE_CCTK_PARAMETERS;

  // Loop over the entire grid (0 to n-1 cells in each direction)
  constexpr auto DI = PointDesc::DI;
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        metric g;
        g.gxx = 0.0;
        g.gxy = 0.0;
        g.gxz = 0.0;
        g.gyy = 0.0;
        g.gyz = 0.0;
        g.gzz = 0.0;

        for (int dk = 0; dk < 2; ++dk) {
          for (int dj = 0; dj < 2; ++dj) {
            for (int di = 0; di < 2; ++di) {
              g.gxx += gxx(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g.gxy += gxy(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g.gxz += gxz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g.gyy += gyy(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g.gyz += gyz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
              g.gzz += gzz(p.I + DI[0] * di + DI[1] * dj + DI[2] * dk);
            }
          }
        }

        g.gxx *= 0.125;
        g.gxy *= 0.125;
        g.gxz *= 0.125;
        g.gyy *= 0.125;
        g.gyz *= 0.125;
        g.gzz *= 0.125;

        prim pv;
        pv.rho = rho(p.I);
        pv.velx = velx(p.I);
        pv.vely = vely(p.I);
        pv.velz = velz(p.I);
        pv.eps = eps(p.I);
        pv.press = press(p.I);

        cons cv;
        prim2con(g, pv, cv);

        dens(p.I) = cv.dens;
        momx(p.I) = cv.momx;
        momy(p.I) = cv.momy;
        momz(p.I) = cv.momz;
        tau(p.I) = cv.tau;
      });
}

} // namespace AsterX
