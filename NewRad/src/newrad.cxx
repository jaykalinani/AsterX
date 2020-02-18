#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <cmath>

namespace NewRad {
using namespace Loop;
using namespace std;

namespace {
template <typename T> constexpr T pow2(const T x) { return x * x; }
} // namespace

// Adapted from BSSN_MoL's files NewRad.F and newrad.h. This code was probably
// originally written by Miguel Alcubierre.
void newrad(const cGH *restrict const cctkGH,
            const CCTK_REAL *restrict const var, CCTK_REAL *restrict const rhs,
            const CCTK_REAL var0, //!< value at infinity
            const CCTK_REAL v0    //!< propagation speed
) {
  DECLARE_CCTK_ARGUMENTS;

  constexpr vect<int, dim> DI{1, 0, 0};
  constexpr vect<int, dim> DJ{0, 1, 0};
  constexpr vect<int, dim> DK{0, 0, 1};

  const CCTK_REAL dx = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz = CCTK_DELTA_SPACE(2);

  const GF3D<const CCTK_REAL, 0, 0, 0> var_(cctkGH, var);
  const GF3D<CCTK_REAL, 0, 0, 0> rhs_(cctkGH, rhs);

  const auto derivx{
      [&](const GF3D<const CCTK_REAL, 0, 0, 0> &u_, const PointDesc &p) {
        const auto I = p.I;
        if (p.NI[0] == 0)
          // interior
          return (u_(I + DI) - u_(I - DI)) / (2 * dx);
        if (p.NI[0] == +1)
          // upper boundary
          return +(3 * u_(I) - 4 * u_(I - DI) + u_(I - 2 * DI)) / (2 * dx);
        if (p.NI[0] == -1)
          // lower boundary
          return -(3 * u_(I) - 4 * u_(I + DI) + u_(I + 2 * DI)) / (2 * dx);
        assert(0);
      }};
  const auto derivy{
      [&](const GF3D<const CCTK_REAL, 0, 0, 0> &u_, const PointDesc &p) {
        const auto I = p.I;
        if (p.NI[1] == 0)
          return (u_(I + DJ) - u_(I - DJ)) / (2 * dy);
        if (p.NI[1] == +1)
          return +(3 * u_(I) - 4 * u_(I - DJ) + u_(I - 2 * DJ)) / (2 * dy);
        if (p.NI[1] == -1)
          return -(3 * u_(I) - 4 * u_(I + DJ) + u_(I + 2 * DJ)) / (2 * dy);
        assert(0);
      }};
  const auto derivz{
      [&](const GF3D<const CCTK_REAL, 0, 0, 0> &u_, const PointDesc &p) {
        const auto I = p.I;
        if (p.NI[2] == 0)
          return (u_(I + DK) - u_(I - DK)) / (2 * dz);
        if (p.NI[2] == +1)
          return +(3 * u_(I) - 4 * u_(I - DK) + u_(I - 2 * DK)) / (2 * dz);
        if (p.NI[2] == -1)
          return -(3 * u_(I) - 4 * u_(I + DK) + u_(I + 2 * DK)) / (2 * dz);
        assert(0);
      }};

  loop_bnd<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // The main part of the boundary condition assumes that we have an
    // outgoing radial wave with some speed v0:
    //
    //    var  =  var0 + u(r-v0*t)/r
    //
    // This implies the following differential equation:
    //
    //    d_t var  =  - v^i d_i var  -  v0 (var - var0) / r
    //
    // where  vi = v0 xi/r

    const CCTK_REAL r = sqrt(pow2(p.x) + pow2(p.y) + pow2(p.z));

    // Find local wave speeds
    const CCTK_REAL vx = v0 * p.x / r;
    const CCTK_REAL vy = v0 * p.y / r;
    const CCTK_REAL vz = v0 * p.z / r;
    // CCTK_REAL const vr = sqrt(pow2(vx) + pow2(vy) + pow2(vz));

    // Derivatives
    const CCTK_REAL varx = derivx(var_, p);
    const CCTK_REAL vary = derivy(var_, p);
    const CCTK_REAL varz = derivz(var_, p);

    // Calculate source term
    rhs_(p.I) =
        -vx * varx - vy * vary - vz * varz - v0 * (var_(p.I) - var0) / r;
  });
}

extern "C" void NewRad_Apply(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_NewRad_Apply;

  newrad(cctkGH, aDD00GF, aDD00_rhsGF, 0, 1);
  newrad(cctkGH, aDD01GF, aDD01_rhsGF, 0, 1);
  newrad(cctkGH, aDD02GF, aDD02_rhsGF, 0, 1);
  newrad(cctkGH, aDD11GF, aDD11_rhsGF, 0, 1);
  newrad(cctkGH, aDD12GF, aDD12_rhsGF, 0, 1);
  newrad(cctkGH, aDD22GF, aDD22_rhsGF, 0, 1);
  newrad(cctkGH, alphaGF, alpha_rhsGF, 1, 1);
  newrad(cctkGH, betU0GF, betU0_rhsGF, 0, 1);
  newrad(cctkGH, betU1GF, betU1_rhsGF, 0, 1);
  newrad(cctkGH, betU2GF, betU2_rhsGF, 0, 1);
  newrad(cctkGH, cfGF, cf_rhsGF, 1, 1);
  newrad(cctkGH, hDD00GF, hDD00_rhsGF, 0, 1);
  newrad(cctkGH, hDD01GF, hDD01_rhsGF, 0, 1);
  newrad(cctkGH, hDD02GF, hDD02_rhsGF, 0, 1);
  newrad(cctkGH, hDD11GF, hDD11_rhsGF, 0, 1);
  newrad(cctkGH, hDD12GF, hDD12_rhsGF, 0, 1);
  newrad(cctkGH, hDD22GF, hDD22_rhsGF, 0, 1);
  newrad(cctkGH, lambdaU0GF, lambdaU0_rhsGF, 0, 1);
  newrad(cctkGH, lambdaU1GF, lambdaU1_rhsGF, 0, 1);
  newrad(cctkGH, lambdaU2GF, lambdaU2_rhsGF, 0, 1);
  newrad(cctkGH, trKGF, trK_rhsGF, 0, 1);
  newrad(cctkGH, vetU0GF, vetU0_rhsGF, 0, 1);
  newrad(cctkGH, vetU1GF, vetU1_rhsGF, 0, 1);
  newrad(cctkGH, vetU2GF, vetU2_rhsGF, 0, 1);
}

} // namespace NewRad
