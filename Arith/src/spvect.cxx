#include "dual.hxx"
#include "spvect.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>

#include <cassert>
#include <sstream>

namespace Arith {

extern "C" void TestSpvect(CCTK_ARGUMENTS) {
  CCTK_INFO("TestSPvect");

  using R = CCTK_REAL;
  using VR = spvect<int, R>;

  const VR n;
  VR x, y, z;
  x.emplace_back(1, 1);
  x.emplace_back(2, 2);
  x.emplace_back(3, 3);
  y.emplace_back(2, 12);
  y.emplace_back(4, 14);
  z.emplace_back(1, 20);
  z.emplace_back(4, 24);

  std::cout << "x=" << x << "\n";
  std::cout << "y=" << y << "\n";
  std::cout << "z=" << z << "\n";

  const R a = -2;
  const R b = -3;

  assert(n == n);

  assert(x == x);
  assert(x + n == x);
  assert(n + x == x);
  assert(+x == n + x);

  assert(x + (-x) == n);
  assert(-(-x) == x);
  assert(x - y == x + (-y));

  assert(y == y);
  assert(x + y == y + x);

  assert(z == z);
  assert((x + y) + z == x + (y + z));

  assert(0 * x == n);
  assert(1 * x == x);
  assert((-1) * x == -x);

  assert(a * x == x * a);
  assert((a + b) * x == a * x + b * x);
  assert((a * b) * x == a * (b * x));
  assert(a * (x + y) == a * x + a * y);

  dual<CCTK_REAL, spvect<int, CCTK_REAL> > dsp;
}
} // namespace Arith
