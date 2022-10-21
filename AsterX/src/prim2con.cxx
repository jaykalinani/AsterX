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
        const smat<CCTK_REAL, 3> g{calc_avg_v2c(gxx, p), calc_avg_v2c(gxy, p),
                                   calc_avg_v2c(gxz, p), calc_avg_v2c(gyy, p),
                                   calc_avg_v2c(gyz, p), calc_avg_v2c(gzz, p)};

        // Interpolate lapse and shift from vertice to center
        const CCTK_REAL lapse = calc_avg_v2c(alp, p);
        const vec<CCTK_REAL, 3> shift{calc_avg_v2c(betax, p),
                                      calc_avg_v2c(betay, p),
                                      calc_avg_v2c(betaz, p)};

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
        prim2con(g, lapse, shift, pv, cv);

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
      });

  /* Initilaize Psi to 0.0 */
  grid.loop_all_device<0, 0, 0>(
      grid.nghostzones, [=] CCTK_DEVICE(const PointDesc &p)
                            CCTK_ATTRIBUTE_ALWAYS_INLINE { Psi(p.I) = 0.0; });
}

} // namespace AsterX
