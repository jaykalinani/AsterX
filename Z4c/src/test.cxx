#warning "TODO"
#include <iostream>

#include "physics.hxx"
#include "tensor.hxx"

#include <cctk.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <random>

namespace Z4c {
using namespace std;

// TODO: Use GoogleTest instead of assert

extern "C" void Z4c_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  mt19937 engine(42);
  uniform_int_distribution<int> dist(-10, 10);
  const auto rand10{[&]() { return double(dist(engine)); }};
  const auto randmat10{[&]() {
    array<array<double, 3>, 3> arr;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        arr[a][b] = rand10();
    return mat3<double>(
        [&](int a, int b) { return arr[min(a, b)][max(a, b)]; });
  }};

  const mat3<double> Z([&](int a, int b) { return double(0); });
  const mat3<double> I([&](int a, int b) { return double(a == b); });
  assert(I != Z);

  for (int n = 0; n < 100; ++n) {
    const mat3<double> A = randmat10();
    const mat3<double> B = randmat10();
    const mat3<double> C = randmat10();
    const double a = rand10();
    const double b = rand10();

    assert((A + B) + C == A + (B + C));
    assert(Z + A == A);
    assert(A + Z == A);
    assert(A + (-A) == Z);
    assert((-A) + A == Z);
    assert(A - B == A + (-B));
    assert(A + B == B + A);

    assert(1 * A == A);
    assert(0 * A == Z);
    assert(-1 * A == -A);
    // assert(mul(a * A, B) == a * mul(A, B));
    assert((a * b) * A == a * (b * A));
    assert(a * (A + B) == a * A + a * B);
    assert((a + b) * A == a * A + b * A);

    // assert(mul(mul(A, B), C) == mul(A, mul(B, C)));
    assert(mul(I, A) == A);
    assert(mul(A, I) == A);
    assert(mul(Z, A) == Z);
    assert(mul(A, Z) == Z);

    assert(Z.det() == 0);
    assert(I.det() == 1);
    assert((a * A).det() == pow3(a) * A.det());

    assert(Z.inv(1) == Z);
    assert(I.inv(1) == I);

    assert(mul(A.inv(1), A) == A.det() * I);
    assert(mul(A, A.inv(1)) == A.det() * I);

    assert((a * A).inv(1) == pow2(a) * A.inv(1));
  }
}

} // namespace Z4c
