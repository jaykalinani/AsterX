#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <cmath>

namespace StaticTrumpet {
using namespace Loop;
using namespace std;

template <typename T> T pow2(const T x) { return x * x; }
template <typename T> T pow4(const T x) {
  const auto x2 = pow2(x);
  return x2 * x2;
}

extern "C" void StaticTrumpet_Initial(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const bool set_metric = CCTK_EQUALS(initial_data, "Static Trumpet");
  const bool set_lapse = CCTK_EQUALS(initial_lapse, "Static Trumpet");
  const bool set_shift = CCTK_EQUALS(initial_shift, "Static Trumpet");
  const bool set_dtshift = CCTK_EQUALS(initial_dtshift, "Static Trumpet");

  const GF3D<CCTK_REAL, 0, 0, 0> gxx_(cctkGH, gxx);
  const GF3D<CCTK_REAL, 0, 0, 0> gxy_(cctkGH, gxy);
  const GF3D<CCTK_REAL, 0, 0, 0> gxz_(cctkGH, gxz);
  const GF3D<CCTK_REAL, 0, 0, 0> gyy_(cctkGH, gyy);
  const GF3D<CCTK_REAL, 0, 0, 0> gyz_(cctkGH, gyz);
  const GF3D<CCTK_REAL, 0, 0, 0> gzz_(cctkGH, gzz);

  const GF3D<CCTK_REAL, 0, 0, 0> kxx_(cctkGH, kxx);
  const GF3D<CCTK_REAL, 0, 0, 0> kxy_(cctkGH, kxy);
  const GF3D<CCTK_REAL, 0, 0, 0> kxz_(cctkGH, kxz);
  const GF3D<CCTK_REAL, 0, 0, 0> kyy_(cctkGH, kyy);
  const GF3D<CCTK_REAL, 0, 0, 0> kyz_(cctkGH, kyz);
  const GF3D<CCTK_REAL, 0, 0, 0> kzz_(cctkGH, kzz);

  const GF3D<CCTK_REAL, 0, 0, 0> alp_(cctkGH, alp);

  // const GF3D<CCTK_REAL, 0, 0, 0> dtalp_(cctkGH, dtalp);

  const GF3D<CCTK_REAL, 0, 0, 0> betax_(cctkGH, betax);
  const GF3D<CCTK_REAL, 0, 0, 0> betay_(cctkGH, betay);
  const GF3D<CCTK_REAL, 0, 0, 0> betaz_(cctkGH, betaz);

  const GF3D<CCTK_REAL, 0, 0, 0> dtbetax_(cctkGH, dtbetax);
  const GF3D<CCTK_REAL, 0, 0, 0> dtbetay_(cctkGH, dtbetay);
  const GF3D<CCTK_REAL, 0, 0, 0> dtbetaz_(cctkGH, dtbetaz);

  constexpr CCTK_REAL rmin = 1.0e-2;
  constexpr CCTK_REAL M = 1;

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    const auto r = sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z));
    const auto r1 = 1 / fmax(rmin, r);

    const auto rho = sqrt(pow2(p.x) * pow2(p.y));
    const auto rho1 = 1 / fmax(rmin, rho);

    // x = r sin theta cos phi
    // y = r sin theta sin phi
    // z = r cos theta

    // rho = r sin theta = sqrt (x^2 + y^2)

    // dx/dr = x / r
    // dy/dr = y / r
    // dz/dr = z / r

    // dx/dtheta =   r cos theta cos phi = x * z / rho
    // dy/dtheta =   r cos theta sin phi = y * z / rho
    // dz/dtheta = - r sin theta         = - rho

    // dx/dphi = - r sin theta sin phi = - y
    // dy/dphi =   r sin theta cos phi =   x
    // dz/dphi =   0

    const auto ex_r = p.x * r1;
    const auto ey_r = p.y * r1;
    const auto ez_r = p.z * r1;
    const auto ex_theta = +p.x * p.z * rho1;
    const auto ey_theta = +p.y * p.z * rho1;
    const auto ez_theta = -rho;
    const auto ex_phi = -p.y;
    const auto ey_phi = +p.x;
    const auto ez_phi = -rho;

    // first index i (cartesian), second index a (spherical)
    const array<array<CCTK_REAL, 3>, 3> dX_dR{{
        {{ex_r, ex_theta, ex_phi}},
        {{ey_r, ey_theta, ey_phi}},
        {{ez_r, ez_theta, ez_phi}},
    }};

    const auto phi = sqrt(1 + M * r1);
    const auto phi4 = pow4(phi);

    // array<array<CCTK_REAL, 3>, 3>
    //   gij{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}};
    // gij[0][0] = phi4;
    // gij[1][1] = phi4;
    // gij[2][2] = phi4;
    // array<array<CCTK_REAL, 3>, 3>
    //   guij{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}};
    // guij[0][0] = 1 / phi4;
    // guij[1][1] = 1 / phi4;
    // guij[2][2] = 1 / phi4;

    // array<array<CCTK_REAL, 3>, 3>
    //   gab{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}};
    // gab[0][0] = phi4;
    // gab[1][1] = pow2(r);
    // gab[2][2] = pow2(rho);
    // array<array<CCTK_REAL, 3>, 3>
    //   guab{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}};
    // guab[0][0] = 1 / phi4;
    // guab[1][1] = pow2(r1);
    // guab[2][2] = pow2(rho1);

    if (set_metric) {
      gxx_(p.I) = phi4;
      gxy_(p.I) = 0;
      gxz_(p.I) = 0;
      gyy_(p.I) = phi4;
      gyz_(p.I) = 0;
      gzz_(p.I) = phi4;
    }

    // Set K_ab in spherical coordinates
    // array<array<CCTK_REAL, 3>, 3>
    //   kab{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}};
    // kab[0][0] = -M * pow2(r1);
    // kab[1][1] = M;
    // kab[2][2] = M * pow2(rho * r1);
    // Raise indices
    array<array<CCTK_REAL, 3>, 3> kuab{{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}};
    kuab[0][0] = -M * pow2(r1) / pow2(phi4);
    kuab[1][1] = M * pow2(r1);
    kuab[2][2] = M * pow2(r1); // M * pow2(rho * r1) * pow2(rho1)
    // Convert to cartesian coordinates
    array<array<CCTK_REAL, 3>, 3> kuij;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        CCTK_REAL t = 0;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b)
            t += dX_dR[i][a] * dX_dR[j][b] * kuab[a][b];
        kuij[i][j] = t;
      }
    // Symmetrize
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < i; ++j) {
        const auto t = (kuij[i][j] + kuij[j][i]) / 2;
        kuij[i][j] = kuij[j][i] = t;
      }
    // Lower indices
    array<array<CCTK_REAL, 3>, 3> kij;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        kij[i][j] = pow2(phi4) * kuij[i][j];

    if (set_metric) {
      kxx_(p.I) = kij[0][0];
      kxy_(p.I) = kij[0][1];
      kxz_(p.I) = kij[0][2];
      kyy_(p.I) = kij[1][1];
      kyz_(p.I) = kij[1][2];
      kzz_(p.I) = kij[2][2];
    }

    if (set_lapse)
      alp_(p.I) = r / (M + r1);

    // dtalp_(p.I) = 0;

    if (set_shift) {
      const auto betar = M * r / pow2(M + r1);
      betax_(p.I) = betar * ex_r;
      betay_(p.I) = betar * ey_r;
      betaz_(p.I) = betar * ez_r;
    }

    if (set_dtshift) {
      dtbetax_(p.I) = 0;
      dtbetay_(p.I) = 0;
      dtbetaz_(p.I) = 0;
    }
  });
}

} // namespace StaticTrumpet
