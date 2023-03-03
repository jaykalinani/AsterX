#include "forms.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <gtest/gtest.h>

#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

namespace TestForms {
using namespace Forms;
using namespace std;

template <typename T> class arbitrary {

  vector<T> vals;
  size_t pos;

  random_device rd;
  default_random_engine gen;
  uniform_int_distribution<> dist;

public:
  arbitrary() : pos(0), gen(rd()), dist(-1000 * 1000, 1000 * 1000) {
    for (int i = 0; i <= 10; ++i)
      vals.push_back(i);
    {
      int n = vals.size();
      for (int i = 0; i < n; ++i)
        vals.push_back(-vals.at(i));
    }
    {
      int n = vals.size();
      for (int i = 0; i < n; ++i)
        vals.push_back(1000 * vals.at(i));
    }
  }

  T next() {
    if (pos < vals.size())
      return vals.at(pos++);
    return dist(gen);
  }
};

template <typename T, int D, int R> void TestVector() {
  arbitrary<T> arb;

  const int n = 100;
  for (int i = 0; i < n; ++i) {

    constexpr int N = form<T, D, R>::N;
    form<T, D, R> x, y, z;
    for (int c = 0; c < N; ++c) {
      x.comps[c] = arb.next();
      y.comps[c] = arb.next();
      z.comps[c] = arb.next();
    }
    const auto zero = form<T, D, R>();

    const T a = arb.next();
    const T b = arb.next();
    const T szero = 0;
    const T sone = 1;

    // Addition:

    // Associativity
    EXPECT_EQ((x + y) + z, x + (y + z));

    // Left zero
    EXPECT_EQ(zero + x, x);
    // Right zero
    EXPECT_EQ(x + zero, x);

    // Commutativity
    EXPECT_EQ(x + y, y + x);

    // Left inverse
    EXPECT_EQ((-x) + x, zero);
    // Right inverse
    EXPECT_EQ(x + (-x), zero);

    // Subtraction
    EXPECT_EQ(x - y, x + (-y));

    // Scaling:

    // Associativity
    EXPECT_EQ((a * b) * x, a * (b * x));

    // Left Unit
    EXPECT_EQ(sone * x, x);
    // Right Unit
    EXPECT_EQ(x * sone, x);

    // Left zero
    EXPECT_EQ(szero * x, zero);
    // Right zero
    EXPECT_EQ(x * szero, zero);

    // Commutativity
    EXPECT_EQ(a * x, x * a);

    if (false) {
      if (a != szero) {
        // Left inverse
        EXPECT_EQ((sone / a) * (a * x), x);
        // Right inverse
        EXPECT_EQ((x * a) * (sone / a), x);

        // Division
        EXPECT_EQ(x / a, x * (sone / a));
      }
    }

    // Vector distributivity
    EXPECT_EQ(a * (x + y), a * x + a * y);

    // Scalar distributivity
    EXPECT_EQ((a + b) * x, a * x + b * x);
  }
}

template <typename T, int D, int Rx, int Ry> void TestWedge() {
  constexpr int R = Rx + Ry;

  arbitrary<T> arb;

  const int n = 100;
  for (int i = 0; i < n; ++i) {

    constexpr int Nx = form<T, D, Rx>::N;
    constexpr int Ny = form<T, D, Ry>::N;
    form<T, D, Rx> x;
    for (int c = 0; c < Nx; ++c)
      x.comps[c] = arb.next();
    form<T, D, Ry> y, y2;
    for (int c = 0; c < Ny; ++c)
      y.comps[c] = arb.next();
    for (int c = 0; c < Ny; ++c)
      y2.comps[c] = arb.next();

    const auto xzero = form<T, D, Rx>();
    const auto yzero = form<T, D, Ry>();
    const auto rzero = form<T, D, R>();

    const T a = arb.next();
    const T b = arb.next();
    const T szero = 0;
    const T sone = 1;

    // Scaling:

    // Associativity
    EXPECT_EQ(a * wedge(x, y), wedge(x, a * y));

    // Left zero
    EXPECT_EQ(wedge(xzero, y), rzero);
    // Right zero
    EXPECT_EQ(wedge(x, yzero), rzero);

    // Distributivity
    EXPECT_EQ(wedge(x, y + y2), wedge(x, y) + wedge(x, y2));

    // Anti-commutativity
    constexpr bool s = (Rx & Ry) & 1;
    EXPECT_EQ(wedge(x, y), bitsign(s) * wedge(y, x));
  }
}

template <typename T, int D, int Rx, int Ry, int Rz> void TestWedge3() {
  constexpr int R = Rx + Ry + Rz;

  arbitrary<T> arb;

  const int n = 100;
  for (int i = 0; i < n; ++i) {

    constexpr int Nx = form<T, D, Rx>::N;
    constexpr int Ny = form<T, D, Ry>::N;
    constexpr int Nz = form<T, D, Rz>::N;
    form<T, D, Rx> x;
    for (int c = 0; c < Nx; ++c)
      x.comps[c] = arb.next();
    form<T, D, Ry> y;
    for (int c = 0; c < Ny; ++c)
      y.comps[c] = arb.next();
    form<T, D, Rz> z;
    for (int c = 0; c < Nz; ++c)
      z.comps[c] = arb.next();

    // Associativity
    EXPECT_EQ(wedge(wedge(x, y), z), wedge(x, wedge(y, z)));
  }
}

extern "C" void TestForms_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_VINFO("Testing...");

  typedef int64_t val_t;

  TestVector<val_t, 0, 0>();
  TestVector<val_t, 1, 0>();
  TestVector<val_t, 1, 1>();
  TestVector<val_t, 2, 0>();
  TestVector<val_t, 2, 1>();
  TestVector<val_t, 2, 2>();
  TestVector<val_t, 3, 0>();
  TestVector<val_t, 3, 1>();
  TestVector<val_t, 3, 2>();
  TestVector<val_t, 3, 3>();

  TestWedge<val_t, 0, 0, 0>();
  TestWedge3<val_t, 0, 0, 0, 0>();

  TestWedge<val_t, 1, 0, 0>();
  TestWedge<val_t, 1, 0, 1>();
  TestWedge<val_t, 1, 1, 0>();
  TestWedge3<val_t, 1, 0, 0, 0>();
  TestWedge3<val_t, 1, 0, 0, 1>();
  TestWedge3<val_t, 1, 0, 1, 0>();
  TestWedge3<val_t, 1, 1, 0, 0>();

  TestWedge<val_t, 2, 0, 0>();
  TestWedge<val_t, 2, 0, 1>();
  TestWedge<val_t, 2, 0, 2>();
  TestWedge<val_t, 2, 1, 0>();
  TestWedge<val_t, 2, 1, 1>();
  TestWedge<val_t, 2, 2, 0>();
  TestWedge3<val_t, 2, 0, 0, 0>();
  TestWedge3<val_t, 2, 0, 0, 1>();
  TestWedge3<val_t, 2, 0, 0, 2>();
  TestWedge3<val_t, 2, 0, 1, 0>();
  TestWedge3<val_t, 2, 0, 1, 1>();
  TestWedge3<val_t, 2, 0, 2, 0>();
  TestWedge3<val_t, 2, 1, 0, 0>();
  TestWedge3<val_t, 2, 1, 0, 1>();
  TestWedge3<val_t, 2, 1, 1, 0>();
  TestWedge3<val_t, 2, 2, 0, 0>();

  TestWedge<val_t, 3, 0, 0>();
  TestWedge<val_t, 3, 0, 1>();
  TestWedge<val_t, 3, 0, 2>();
  TestWedge<val_t, 3, 0, 3>();
  TestWedge<val_t, 3, 1, 0>();
  TestWedge<val_t, 3, 1, 1>();
  TestWedge<val_t, 3, 1, 2>();
  TestWedge<val_t, 3, 2, 0>();
  TestWedge<val_t, 3, 2, 1>();
  TestWedge<val_t, 3, 3, 0>();
  TestWedge3<val_t, 3, 0, 0, 0>();
  TestWedge3<val_t, 3, 0, 0, 1>();
  TestWedge3<val_t, 3, 0, 0, 2>();
  TestWedge3<val_t, 3, 0, 0, 3>();
  TestWedge3<val_t, 3, 0, 1, 0>();
  TestWedge3<val_t, 3, 0, 1, 1>();
  TestWedge3<val_t, 3, 0, 1, 2>();
  TestWedge3<val_t, 3, 0, 2, 0>();
  TestWedge3<val_t, 3, 0, 2, 1>();
  TestWedge3<val_t, 3, 1, 0, 0>();
  TestWedge3<val_t, 3, 1, 0, 1>();
  TestWedge3<val_t, 3, 1, 0, 2>();
  TestWedge3<val_t, 3, 1, 1, 0>();
  TestWedge3<val_t, 3, 1, 1, 1>();
  TestWedge3<val_t, 3, 2, 0, 0>();
  TestWedge3<val_t, 3, 2, 0, 1>();
  TestWedge3<val_t, 3, 2, 1, 0>();
  TestWedge3<val_t, 3, 3, 0, 0>();

  CCTK_VINFO("Done.");
}

} // namespace TestForms
