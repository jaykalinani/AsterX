#include "simd.hxx"

#include <fixmath.hxx> // include this before <cctk.h>
#include <cctk.h>

#include <cassert>

namespace Arith {

void TestSIMD() {
  // nvcc V11.1.74 doesn't accept this as "constexpr" values
#ifndef __CUDACC__
  typedef simd<CCTK_REAL> realv;

  realv x;
  realv y = 0;
  realv z = zero<realv>();

  assert(all(y == 0));
  assert(all(y == z));
#endif
}

} // namespace Arith
