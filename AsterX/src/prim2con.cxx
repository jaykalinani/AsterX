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
  grid.loop_int_device<1, 1, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        // Interpolate metric terms from vertices to center
        metric g;
        g.gxx = calc_avg_v2c(gxx, p);
        g.gxy = calc_avg_v2c(gxy, p);
        g.gxz = calc_avg_v2c(gxz, p);
        g.gyy = calc_avg_v2c(gyy, p);
        g.gyz = calc_avg_v2c(gyz, p);
        g.gzz = calc_avg_v2c(gzz, p);
        
	// Interpolate lapse and shift from vertice to center
	lapse l;
        l.alp = calc_avg_v2c(alp, p);
	shift shft;
        shft.betax = calc_avg_v2c(betax, p);
        shft.betay = calc_avg_v2c(betay, p);
        shft.betaz = calc_avg_v2c(betaz, p);

        prim pv;
        pv.rho = rho(p.I);
        pv.velx = velx(p.I);
        pv.vely = vely(p.I);
        pv.velz = velz(p.I);
        pv.eps = eps(p.I);
        pv.press = press(p.I);
        pv.Bvecx = Bvecx(p.I);
        pv.Bvecy = Bvecy(p.I);
        pv.Bvecz = Bvecz(p.I);

        cons cv;
        prim2con(g, l, shft, pv, cv);

        dens(p.I) = cv.dens;
        momx(p.I) = cv.momx;
        momy(p.I) = cv.momy;
        momz(p.I) = cv.momz;
        tau(p.I) = cv.tau;
        dBx(p.I) = cv.dBvecx;
        dBy(p.I) = cv.dBvecy;
        dBz(p.I) = cv.dBvecz;

        saved_rho(p.I) = pv.rho;
        saved_velx(p.I) = pv.velx;
        saved_vely(p.I) = pv.vely;
        saved_velz(p.I) = pv.velz;
        saved_eps(p.I) = pv.eps;
        saved_Bvecx(p.I) = pv.Bvecx;
        saved_Bvecy(p.I) = pv.Bvecy;
        saved_Bvecz(p.I) = pv.Bvecz;
      });

  /* Initialize Avec_x, Ex */
  grid.loop_int_device<1, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Ex(p.I) = 0.0;
        Avec_x(p.I) = 0.0; // TODO fix this
      });
  /* Initialize Avec_y, Ey */
  grid.loop_int_device<0, 1, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Ey(p.I) = 0.0;
        Avec_y(p.I) = 0.0; // TODO fix this
      });
  /* Initialize Avec_z, Ez */
  grid.loop_int_device<0, 0, 1>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Ez(p.I) = 0.0;
        Avec_z(p.I) = 0.0; // TODO fix this
      });
  /* Initilaize Psi */
  grid.loop_int_device<0, 0, 0>(
      grid.nghostzones,
      [=] CCTK_DEVICE(const PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
        Psi(p.I) = 0.0; // TODO fix this
      });
}

} // namespace AsterX
