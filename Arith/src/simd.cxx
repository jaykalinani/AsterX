#include "simd.hxx"

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

  realv a = 2;
  realv b = 3;
  realv c = 4;

  realv r0 = muladd(a, b, c);
  realv r1 = mulsub(a, b, c);
  realv r2 = negmuladd(a, b, c);
  realv r3 = negmulsub(a, b, c);
  assert(all(r0 == muladd(2, 3, 4)));
  assert(all(r1 == mulsub(2, 3, 4)));
  assert(all(r2 == negmuladd(2, 3, 4)));
  assert(all(r3 == negmulsub(2, 3, 4)));
#endif
}

} // namespace Arith
