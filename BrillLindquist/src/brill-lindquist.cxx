#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace BrillLindquist {
using namespace Loop;

template <typename T> constexpr T pow2(const T x) { return x * x; }
template <typename T> constexpr T pow4(const T x) { return pow2(x) * pow2(x); }

template <typename T> constexpr T expand(const T x) {
  constexpr T x0 = 5;
  return x;
  // return sinh(x / x0) * x0;
}

template <typename T> constexpr T smooth(const T r) {
  DECLARE_CCTK_PARAMETERS;
  const T rmin = min_radius;
  return fmax(rmin, r);
}

template <typename T> constexpr T psi(const T x, const T y, const T z) {
  DECLARE_CCTK_PARAMETERS;
  const T m = mass;
  const T r = smooth(
      sqrt(pow2(expand(x - x0)) + pow2(expand(y - y0)) + pow2(expand(z - z0))));
  const T phi = 1 + m / (2 * r);
  return pow4(phi);
}

template <typename T> constexpr T lapse(const T x, const T y, const T z) {
  DECLARE_CCTK_PARAMETERS;
  const T m = mass;
  const T r = smooth(
      sqrt(pow2(expand(x - x0)) + pow2(expand(y - y0)) + pow2(expand(z - z0))));
  const T lapse1 = (1 - m / (2 * r)) / (1 + m / (2 * r));
  return (1 + lapse1) / 2; // average
}

constexpr array<int, 3> CarpetX_widths{17, 5, 5};

bool CarpetX(int i, int j, int k) {
  // Logo originally generated via
  // <http://patorjk.com/software/taag/#p=display&f=Alphabet&t=CAR%20PET%20X>,
  // then hand-modified
  const vector<string> CAR{
      // width: 5 + 5 + 5, height: 5, XZ
      " CCCC  AAA  RRRR ", //
      "C     A   A R   R", //
      "C     AAAAA RRRR ", //
      "C     A   A R R  ", //
      " CCCC A   A R  RR", //
  };
  const vector<string> PET{
      // widht: 5 + 5 + 5, height: 5, XY
      "PPPP  EEEEE TTTTT", //
      "P   P E       T  ", //
      "PPPP  EEEE    T  ", //
      "P     E       T  ", //
      "P     EEEEE   T  ", //
  };
  const vector<string> X{
      // width: 5, height: 5, YZ
      "X   X ", //
      " X X  ", //
      "  X   ", //
      " X X  ", //
      "X   X ", //
  };
  assert(int(CAR.size()) == CarpetX_widths[2]);
  assert(int(PET.size()) == CarpetX_widths[1]);
  assert(int(X.size()) == CarpetX_widths[2]);
  for (const auto &car : CAR)
    assert(int(car.size()) == CarpetX_widths[0]);
  for (const auto &pet : PET)
    assert(int(pet.size()) == CarpetX_widths[0]);
  for (const auto &x : X)
    assert(int(x.size()) == CarpetX_widths[1]);
  assert(i >= 0 && i < CarpetX_widths[0]);
  assert(j >= 0 && j < CarpetX_widths[1]);
  assert(k >= 0 && k < CarpetX_widths[2]);
  return CAR.at(k).at(i) != ' ' && PET.at(j).at(i) != ' ' &&
         X.at(k).at(j) != ' ';
}

inline CCTK_REAL psi(const PointDesc &p) { return psi(p.x, p.y, p.z); }
inline CCTK_REAL lapse(const PointDesc &p) { return lapse(p.x, p.y, p.z); }

extern "C" void BrillLindquist_initial_data(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_BrillLindquist_initial_data;
  DECLARE_CCTK_PARAMETERS;

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

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxx_(p.I) = psi(p); });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxy_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gxz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gyy_(p.I) = psi(p); });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gyz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { gzz_(p.I) = psi(p); });

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxx_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxy_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kxz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kyy_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kyz_(p.I) = 0; });
  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) { kzz_(p.I) = 0; });
}

extern "C" void BrillLindquist_initial_lapse(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_BrillLindquist_initial_lapse;
  DECLARE_CCTK_PARAMETERS;

  const GF3D<CCTK_REAL, 0, 0, 0> alp_(cctkGH, alp);

  loop_all<0, 0, 0>(cctkGH, [&](const PointDesc &p) {
    // alp_(p.I) = 1 / psi(p);
    alp_(p.I) = lapse(p);
  });
}

} // namespace BrillLindquist
