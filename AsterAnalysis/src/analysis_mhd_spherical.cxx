#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>
#include <limits>

#include "aster_utils.hxx"

namespace AsterAnalysis {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void AsterAnalysis_MHD_Spherical(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterAnalysis_MHD_Spherical;
  DECLARE_CCTK_PARAMETERS;

  DECLARE_CCTK_ARGUMENTSX_AsterAnalysis_MHD_Spherical;
  DECLARE_CCTK_PARAMETERS;

  /* aliases – metric & B are read‑only; outputs are direct arrays */
  const smat<GF3D2<const CCTK_REAL>, 3> gf_g{gxx, gxy, gxz, gyy, gyz, gzz};
  const vec<GF3D2<const CCTK_REAL>, 3> gf_B{Bvecx, Bvecy, Bvecz};

  /* pre‑computed inverse grid spacings (uniform) */
  const CCTK_REAL det_eps = std::numeric_limits<CCTK_REAL>::epsilon();

  /* identity rotation (toroidalAxis = +z) */
  constexpr smat<CCTK_REAL, 3> R{{{1,0,0},{0,1,0},{0,0,1}}};

  /* main loop over all zone‑centred points */
  grid.loop_all_device<1,1,1>(grid.nghostzones,
    [=] CCTK_DEVICE (const PointDesc& p) CCTK_ATTRIBUTE_ALWAYS_INLINE {

      /* coordinates in (possibly rotated) Cartesian frame */
      const CCTK_REAL newX = x(p.I);      //Spritz: R ⋅ (x - origin), but origin = 0 and R = I
      const CCTK_REAL newY = y(p.I);
      const CCTK_REAL newZ = z(p.I);

      /* metric at centre */
      const smat<CCTK_REAL,3> g([&](int i,int j) ARITH_INLINE {
        return gf_g(i,j)(p.I);
      });
      const CCTK_REAL detg = calc_det(g);
      const CCTK_REAL sqrt_detg = sqrt(detg);
      const smat<CCTK_REAL,3> ug = calc_inv(g, detg);

      /* magnetic field – up */
      const vec<CCTK_REAL,3> B_up{gf_B(0)(p.I), gf_B(1)(p.I), gf_B(2)(p.I)}; // already in lab frame
      /* magnetic field – down (g_{ij} B^j) */
      const vec<CCTK_REAL,3> B_low = calc_contraction(g, B_up);

      /* B^2 & |B| */
      const CCTK_REAL B2big = calc_contraction(B_up, B_low);
      const CCTK_REAL B_norm_local = sqrt(B2big);

      /* spherical radii */
      const CCTK_REAL at_r = sqrt(newX*newX + newY*newY);
      const CCTK_REAL at_rho = sqrt(at_r*at_r + newZ*newZ);

      /* outputs default to zero */
      CCTK_REAL Bphi_up = 0.0, Btheta_up = 0.0, Brho_up = 0.0;
      CCTK_REAL Blowphi = 0.0, Bpol_norm = 0.0, Btor_norm = 0.0;

      if (fabs(at_r) > det_eps && fabs(at_rho) > det_eps) {
        /* covariant components in rotated frame (R = I) */
        const CCTK_REAL Blowx = B_low(0);
        const CCTK_REAL Blowy = B_low(1);
        const CCTK_REAL Blowz = B_low(2);

        /* contravariant components already B_up */
        const CCTK_REAL Bupx = B_up(0);
        const CCTK_REAL Bupy = B_up(1);
        const CCTK_REAL Bupz = B_up(2);

        /* φ‑component */
        Bphi_up = (Bupy * newX - Bupx * newY) / (at_r * at_r);
        Blowphi =  Blowy * newX - Blowx * newY;

        /* ρ‑component */
        Brho_up = (Bupx * newX + Bupy * newY + Bupz * newZ) / at_rho;

        /* θ‑component */
        CCTK_REAL tmp = (Bupx * newX + Bupy * newY) * newZ / at_r - Bupz * at_r;
        Btheta_up = tmp / (at_rho * at_rho);

        /* norms */
        Btor_norm = sqrt(fabs(Blowphi * Bphi_up));
        Bpol_norm = sqrt(fabs(B_norm_local*B_norm_local - Blowphi * Bphi_up));
      }

      /* write to grid functions */
      Bphi(p.I)              = Bphi_up;
      Btheta(p.I)            = Btheta_up;
      Brho(p.I)              = Brho_up;
      B_Poloidal_norm(p.I)   = Bpol_norm;
      B_Toroidal_norm(p.I)   = Btor_norm;

  

}

} // namespace AsterAnalysis
