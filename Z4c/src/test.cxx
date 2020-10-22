#include "derivs.hxx"
#include "physics.hxx"
#include "tensor.hxx"

#include <cctk.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <random>

namespace Z4c {
using namespace std;

// TODO: Use GoogleTest instead of assert

extern "C" void Z4c_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

#ifdef __CUDACC_

  // Test tensors

  mt19937 engine(42);
  uniform_int_distribution<int> dist(-10, 10);
  const auto rand10{[&]() { return double(dist(engine)); }};
  const auto randmat10{[&]() {
    array<array<double, 3>, 3> arr;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
        arr[a][b] = rand10();
    return mat3<double, DN, DN>(
        [&](int a, int b) { return arr[min(a, b)][max(a, b)]; });
  }};

  const mat3<double, DN, DN> Z([&](int a, int b) { return double(0); });
  const mat3<double, DN, DN> I([&](int a, int b) { return double(a == b); });
  assert(I != Z);
  const mat3<double, UP, UP> Zup([&](int a, int b) { return double(0); });
  const mat3<double, UP, UP> Iup([&](int a, int b) { return double(a == b); });
  assert(Iup != Zup);

  for (int n = 0; n < 100; ++n) {
    const mat3<double, DN, DN> A = randmat10();
    const mat3<double, DN, DN> B = randmat10();
    const mat3<double, DN, DN> C = randmat10();
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
    // DNUP  assert(mul(I, A) == A);
    // DNUP  assert(mul(A, I) == A);
    // DNUP  assert(mul(Z, A) == Z);
    // DNUP  assert(mul(A, Z) == Z);

    assert(Z.det() == 0);
    assert(I.det() == 1);
    assert((a * A).det() == pow3(a) * A.det());

    assert(Z.inv(1) == Zup);
    assert(I.inv(1) == Iup);

    // DNUP assert(mul(A.inv(1), A) == A.det() * Iup);
    // DNUP assert(mul(A, A.inv(1)) == A.det() * Iup);

    assert((a * A).inv(1) == pow2(a) * A.inv(1));
  }

  // Test derivatives

  static_assert(deriv_order % 2 == 0, "");
  constexpr int required_ghosts = deriv_order / 2 + 1;
  const double eps = 1.0e-12;

  // deriv
  for (int order = 0; order <= deriv_order; ++order) {
    // CCTK_VINFO("Testing deriv order=%d", order);
    array<double, 2 * required_ghosts + 7> arr;
    for (size_t i = 0; i < arr.size(); ++i)
      arr[i] = NAN;
    double *const var = &arr[arr.size() / 2];
    for (int i = -deriv_order / 2; i <= deriv_order / 2; ++i)
      var[i] = pown(i, order);
    const double expected = order == 1 ? 1 : 0;
    const double found = deriv1d(var, 1, 1.0);
    assert(fabs(found - expected) <= eps);
  }

  // deriv (upwind)
  for (int order = 0; order <= deriv_order; ++order) {
    for (int sign = 0; sign <= 1; ++sign) {
      // CCTK_VINFO("Testing deriv (upwind) order=%d sign=%d", order, sign);
      array<double, 2 * required_ghosts + 7> arr;
      for (size_t i = 0; i < arr.size(); ++i)
        arr[i] = NAN;
      double *const var = &arr[arr.size() / 2];
      if (sign)
        for (int i = -deriv_order / 2 - 1; i <= deriv_order / 2 - 1; ++i)
          var[i] = pown(i, order);
      else
        for (int i = -deriv_order / 2 + 1; i <= deriv_order / 2 + 1; ++i)
          var[i] = pown(i, order);
      const double expected = order == 1 ? 1 : 0;
      const double found = deriv1d_upwind(var, 1, sign, 1.0);
      assert(fabs(found - expected) <= eps);
    }
  }

  // deriv2
  for (int order = 0; order <= deriv_order; ++order) {
    // CCTK_VINFO("Testing deriv2 order=%d", order);
    array<double, 2 * required_ghosts + 7> arr;
    for (size_t i = 0; i < arr.size(); ++i)
      arr[i] = NAN;
    double *const var = &arr[arr.size() / 2];
    for (int i = -deriv_order / 2; i <= deriv_order / 2; ++i)
      var[i] = pown(i, order);
    const double expected = order == 2 ? 2 : 0;
    const double found = deriv2_1d(var, 1, 1.0);
    assert(fabs(found - expected) <= eps);
  }

  // deriv2 (mixed)
  for (int orderj = 0; orderj <= deriv_order; ++orderj) {
    for (int orderi = 0; orderi <= deriv_order; ++orderi) {
      // CCTK_VINFO("Testing deriv2 (mixed) order=%d,%d", orderi, orderj);
      array<array<double, 2 * required_ghosts + 7>, 2 * required_ghosts + 7>
          arr;
      for (size_t j = 0; j < arr.size(); ++j)
        for (size_t i = 0; i < arr[j].size(); ++i)
          arr[j][i] = NAN;
      const int di = 1;
      const int dj = arr[0].size();
      double *const var = &arr[arr.size() / 2][arr[0].size() / 2];
      for (int j = -deriv_order / 2; j <= deriv_order / 2; ++j)
        for (int i = -deriv_order / 2; i <= deriv_order / 2; ++i)
          var[j * dj + i * di] = pown(i, orderi) * pown(j, orderj);
      const double expected = orderi == 1 && orderj == 1 ? 1 : 0;
      const double found = deriv2_2d(var, di, dj, 1.0, 1.0);
      assert(fabs(found - expected) <= eps);
    }
  }

  // deriv (dissipation)
  for (int order = 0; order <= deriv_order + 2; ++order) {
    // CCTK_VINFO("Testing deriv (dissipation) order=%d", order);
    array<double, 2 * required_ghosts + 7> arr;
    for (size_t i = 0; i < arr.size(); ++i)
      arr[i] = NAN;
    double *const var = &arr[arr.size() / 2];
    for (int i = -deriv_order / 2 - 1; i <= deriv_order / 2 + 1; ++i)
      var[i] = pown(i, order);
    const double expected = order == deriv_order + 2 ? factorial(order) : 0;
    const double found = deriv1d_diss(var, 1, 1.0);
    assert(fabs(found - expected) <= eps);
  }

#endif
}

} // namespace Z4c
