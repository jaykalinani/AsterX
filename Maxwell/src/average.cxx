#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

namespace Maxwell {

namespace {
template <typename T, int CI, int CJ, int CK>
T average(const Loop::GF3D<const T, CI, CJ, CK> &var,
          const Loop::PointDesc &p) {
  const auto DI = Loop::vect<int, Loop::dim>::unit(0);
  const auto DJ = Loop::vect<int, Loop::dim>::unit(1);
  const auto DK = Loop::vect<int, Loop::dim>::unit(2);
  T res = 0;
  for (int k = 0; k < 2 - CK; ++k)
    for (int j = 0; j < 2 - CJ; ++j)
      for (int i = 0; i < 2 - CI; ++i)
        res += var(p.I + i * DI + j * DJ + k * DK);
  return res / ((2 - CI) * (2 - CJ) * (2 - CK));
}
} // namespace

extern "C" void Maxwell_Average(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_Maxwell_Average;
  DECLARE_CCTK_PARAMETERS;

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> phi_(cctkGH, phi);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ax_(cctkGH, ax);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ay_(cctkGH, ay);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> az_(cctkGH, az);

  const Loop::GF3D<const CCTK_REAL, 1, 0, 0> ex_(cctkGH, ex);
  const Loop::GF3D<const CCTK_REAL, 0, 1, 0> ey_(cctkGH, ey);
  const Loop::GF3D<const CCTK_REAL, 0, 0, 1> ez_(cctkGH, ez);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> byz_(cctkGH, byz);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> bzx_(cctkGH, bzx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> bxy_(cctkGH, bxy);

  const Loop::GF3D<const CCTK_REAL, 0, 1, 1> curlayz_(cctkGH, curlayz);
  const Loop::GF3D<const CCTK_REAL, 1, 0, 1> curlazx_(cctkGH, curlazx);
  const Loop::GF3D<const CCTK_REAL, 1, 1, 0> curlaxy_(cctkGH, curlaxy);

  const Loop::GF3D<const CCTK_REAL, 0, 0, 0> dive_(cctkGH, dive);

  const Loop::GF3D<const CCTK_REAL, 1, 1, 1> divb_(cctkGH, divb);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgphi_(cctkGH, avgphi);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgax_(cctkGH, avgax);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgay_(cctkGH, avgay);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgaz_(cctkGH, avgaz);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgex_(cctkGH, avgex);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgey_(cctkGH, avgey);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgez_(cctkGH, avgez);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgbyz_(cctkGH, avgbyz);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgbzx_(cctkGH, avgbzx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgbxy_(cctkGH, avgbxy);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgcurlayz_(cctkGH, avgcurlyz);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgcurlazx_(cctkGH, avgcurlzx);
  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgcurlaxy_(cctkGH, avgcurlxy);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgdive_(cctkGH, avgdive);

  const Loop::GF3D<CCTK_REAL, 1, 1, 1> avgdivb_(cctkGH, avgdivb);

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgphi_(p.I) = average(phi_, p);
  });

  Loop::loop_int<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { avgax_(p.I) = average(ax_, p); });
  Loop::loop_int<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { avgay_(p.I) = average(ay_, p); });
  Loop::loop_int<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { avgaz_(p.I) = average(az_, p); });

  Loop::loop_int<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { avgex_(p.I) = average(ex_, p); });
  Loop::loop_int<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { avgey_(p.I) = average(ey_, p); });
  Loop::loop_int<1, 1, 1>(
      cctkGH, [&](const Loop::PointDesc &p) { avgez_(p.I) = average(ez_, p); });

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgbyz_(p.I) = average(byz_, p);
  });
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgbzx_(p.I) = average(bzx_, p);
  });
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgbxy_(p.I) = average(bxy_, p);
  });

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgcurlayz_(p.I) = average(curlayz_, p);
  });
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgcurlazx_(p.I) = average(curlazx_, p);
  });
  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgcurlaxy_(p.I) = average(curlaxy_, p);
  });

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgdive_(p.I) = average(dive_, p);
  });

  Loop::loop_int<1, 1, 1>(cctkGH, [&](const Loop::PointDesc &p) {
    avgdivb_(p.I) = average(divb_, p);
  });
}

} // namespace Maxwell
