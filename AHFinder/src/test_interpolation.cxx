#include <cctk.h>
#include <cctk_Arguments_Checked.h>

#include <array>
#include <cassert>
#include <vector>

namespace AHFinder {
using namespace std;

extern "C" void AHFinder_test_interpolation(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS_AHFinder_test_interpolation;

  const vector<CCTK_INT> all_operations{0, 1, 2, 3, 11, 12, 13, 22, 23, 33};

  const vector<CCTK_INT> all_varinds{
      CCTK_VarIndex("Coordinates::coordx"),
      CCTK_VarIndex("Coordinates::coordy"),
      CCTK_VarIndex("Coordinates::coordz"),
  };

  const int nvars = all_varinds.size() * all_operations.size();

  vector<CCTK_INT> varinds;
  vector<CCTK_INT> operations;
  varinds.reserve(nvars);
  operations.reserve(nvars);
  for (const auto op : all_operations) {
    for (const auto var : all_varinds) {
      varinds.push_back(var);
      operations.push_back(op);
    }
  }
  assert(int(varinds.size()) == nvars);
  assert(int(operations.size()) == nvars);

  constexpr int npoints = 4;
  const array<vector<CCTK_REAL>, 3> coords{
      vector<CCTK_REAL>{0.1, 0.2, 0.3, 0.1},
      vector<CCTK_REAL>{0.1, 0.2, 0.3, 0.2},
      vector<CCTK_REAL>{0.1, 0.2, 0.3, 0.3},
  };
  for (int d = 0; d < 3; ++d)
    assert(coords[d].size() == npoints);

  vector<vector<CCTK_REAL> > results(nvars);
  vector<CCTK_REAL *> resultptrs(nvars);
  for (int n = 0; n < nvars; ++n) {
    results.at(n).resize(npoints);
    resultptrs.at(n) = results.at(n).data();
  }

  Interpolate(cctkGH, npoints, coords[0].data(), coords[1].data(),
              coords[2].data(), nvars, varinds.data(), operations.data(),
              resultptrs.data());

  const auto chop{[](const auto x) { return fabs(x) <= 1.0e-12 ? 0 : x; }};
  for (int v = 0; v < nvars; ++v) {
    for (int n = 0; n < npoints; ++n) {
      const int d = v % all_varinds.size();
      assert(d >= 0 && d < 3);
      const int op = operations.at(v);
      switch (op) {
      case 0:
        assert(chop(results.at(v).at(n) - coords[d].at(n)) == 0);
        break;
      case 1:
      case 2:
      case 3:
        assert(chop(results.at(v).at(n) - (op - 1 == d ? 1 : 0)) == 0);
        break;
      case 11:
      case 12:
      case 13:
      case 22:
      case 23:
      case 33:
        assert(chop(results.at(v).at(n)) == 0);
        break;
      default:
        assert(0);
      }
    }
  }
}

} // namespace AHFinder
