#include <loop.hxx>

#include <cctk.h>
#include <cctk_Arguments_Checked.h>
#include <cctk_Parameters.h>

#include <vector>

namespace AHFinder {
using namespace std;

extern "C" void AHFinder_interpolate(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_interpolate;
  DECLARE_CCTK_PARAMETERS;

  constexpr int npoints = 4;
  const vector<CCTK_REAL> coordsx{0.1, 0.2, 0.3, 0.1};
  const vector<CCTK_REAL> coordsy{0.1, 0.2, 0.3, 0.2};
  const vector<CCTK_REAL> coordsz{0.1, 0.2, 0.3, 0.3};
  constexpr int nvars = 3;
  const vector<CCTK_INT> varinds{
      CCTK_VarIndex("Coordinates::coordx"),
      CCTK_VarIndex("Coordinates::coordy"),
      CCTK_VarIndex("Coordinates::coordz"),
  };
  const vector<CCTK_INT> operations{
      0,
      0,
      0,
  };
  vector<vector<CCTK_REAL> > results(nvars);
  vector<CCTK_REAL *> resultptrs(nvars);
  for (int n = 0; n < nvars; ++n) {
    results.at(n).resize(npoints);
    resultptrs.at(n) = results.at(n).data();
  }
  Interpolate(cctkGH, npoints, coordsx.data(), coordsy.data(), coordsz.data(),
              nvars, varinds.data(), operations.data(), resultptrs.data());

#pragma omp critical
  for (int v = 0; v < nvars; ++v)
    for (int n = 0; n < npoints; ++n)
      CCTK_VINFO("%d[%d]=%g", v, n, double(results.at(v).at(n)));
}

} // namespace AHFinder
