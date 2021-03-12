#include "simt.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

namespace SIMT {
using std::abs;

template <typename T, std::size_t N>
std::ostream &output_array(std::ostream &os, const std::array<T, N> &arr) {
  os << "[";
  for (size_t n = 0; n < arr.size(); ++n) {
    if (n != 0)
      os << ",";
    os << arr[n];
  }
  os << "]";
  return os;
}

std::size_t num_tests = 1234;
std::size_t num_inaccurate_tests = 1234;
std::size_t num_failed_tests = 1234;

#define TEST_EQ(X, Y) test_eq(__FILE__, __LINE__, X, Y)
#define TEST_NE(X, Y) test_ne(__FILE__, __LINE__, X, Y)

template <typename T>
std::enable_if_t<std::is_arithmetic_v<T>, void>
test_eq(const char *file, int line, const T &x, const T &y) {
  ++num_tests;
  const std::equal_to<T> equal;
  if (equal(x, y))
    return;
  const auto delta = maxval(abs(x - y) / max(T(1), max(abs(x), abs(y))));
  const auto eps = 10 * std::numeric_limits<T>::epsilon();
  if (delta <= eps) {
    ++num_inaccurate_tests;
    return;
  }
  std::cout << file << ":" << line
            << ": Test failure in: Expected x==y, found\n"
            << "  x=" << x << "\n"
            << "  y=" << y << "\n"
            << "  Δ=" << (y - x) << "\n"
            << "  |Δ|=" << delta << "\n"
            << "  ε=" << eps << "\n";
}

template <typename T>
void test_eq(const char *file, int line, const simt<T> &x, const simt<T> &y) {
  ++num_tests;
  const std::equal_to<simt<T> > equal;
  if (equal(x, y))
    return;
  const auto delta = maxval(abs(x - y) / max(simt<T>(1), max(abs(x), abs(y))));
  const auto eps = 10 * std::numeric_limits<T>::epsilon();
  if (delta <= eps) {
    ++num_inaccurate_tests;
    return;
  }
  std::cout << file << ":" << line
            << ": Test failure in: Expected x==y, found\n"
            << "  x=" << x << "\n"
            << "  y=" << y << "\n"
            << "  Δ=" << (y - x) << "\n"
            << "  |Δ|=" << delta << "\n"
            << "  ε=" << eps << "\n";
}

template <typename T, std::size_t N>
void test_eq(const char *file, int line, const std::array<T, N> &x,
             const std::array<T, N> &y) {
  ++num_tests;
  const std::equal_to<std::array<T, N> > equal;
  if (equal(x, y))
    return;
  ++num_failed_tests;
  std::cout << file << ":" << line
            << ": Test failure in: Expected x==y, found\n"
            << "  x=";
  output_array(std::cout, x);
  std::cout << "\n"
            << "  y=";
  output_array(std::cout, y);
  std::cout << "\n";
}

template <typename T>
void test_ne(const char *file, int line, const simt<T> &x, const simt<T> &y) {
  ++num_tests;
  const std::equal_to<simt<T> > equal;
  if (!equal(x, y))
    return;
  const auto delta = maxval(abs(x - y) / max(simt<T>(1), max(abs(x), abs(y))));
  const auto eps = 10 * std::numeric_limits<T>::epsilon();
  if (!(delta <= eps)) {
    ++num_inaccurate_tests;
    return;
  }
  ++num_failed_tests;
  std::cout << file << ":" << line
            << ": Test failure in: Expected x!=y, found\n"
            << "  x=" << x << "\n"
            << "  y=" << y << "\n"
            << "  Δ=" << (y - x) << "\n"
            << "  |Δ|=" << delta << "\n"
            << "  ε=" << eps << "\n";
}

template <typename T> void run_tests() {
  std::random_device dev;
  std::default_random_engine gen(dev());
  std::normal_distribution<T> dist(0, 1);

  for (int iter = 0; iter < 100; ++iter) {

    {
      simt<T> w0 CCTK_ATTRIBUTE_UNUSED;
      simt<T> w1 CCTK_ATTRIBUTE_UNUSED{};
      simt<T> w2 CCTK_ATTRIBUTE_UNUSED = {};
    }

    const simt<T> n = 0;
    const simt<T> e = 1;
    const simt<T> x = random(simt<T>{}, gen, dist);
    const simt<T> y = random(simt<T>{}, gen, dist);
    const simt<T> z = random(simt<T>{}, gen, dist);
    const T a = dist(gen);
    const T b = dist(gen);

    TEST_EQ(+x, x);
    TEST_EQ(x + n, x);
    TEST_EQ(n + x, x);
    TEST_EQ(x + y, y + x);
    TEST_EQ((x + y) + z, x + (y + z));

    TEST_EQ(-(-x), x);
    TEST_EQ(n - x, -x);
    TEST_EQ(x - x, n);
    TEST_EQ(-(x + y), (-x) + (-y));

    TEST_EQ(T(0) * x, n);
    TEST_EQ(T(1) * x, x);
    TEST_EQ(a * n, n);
    TEST_EQ(a * x, x * a);
    TEST_EQ((a * b) * x, a * (b * x));
    TEST_EQ(a * (x + y), a * x + a * y);
    TEST_EQ(x * (1 / a), x / a);
    TEST_EQ(x * a * (1 / a), x);
    TEST_EQ(T(-1) * x, -x);

    TEST_EQ(abs(abs(x)), abs(x));
    TEST_EQ(abs(x), abs(-x));
    TEST_EQ(max(n, n), n);
    TEST_EQ(max(y, x), max(x, y));
    TEST_EQ(max(max(x, y), z), max(x, max(y, z)));
    TEST_EQ(max(x, y) + z, max(x + z, y + z));
    TEST_EQ(max(x, -x), abs(x));
    TEST_EQ(min(n, n), n);
    TEST_EQ(min(y, x), min(x, y));
    TEST_EQ(min(min(x, y), z), min(x, min(y, z)));
    TEST_EQ(min(x, y) + z, min(x + z, y + z));
    TEST_EQ(min(x, -x), -abs(x));

    TEST_EQ(abs(a * x), abs(a) * abs(x));
    TEST_EQ(max(abs(a) * x, abs(a) * y), abs(a) * max(x, y));
    TEST_EQ(min(abs(a) * x, abs(a) * y), abs(a) * min(x, y));

    TEST_NE(e, n);
    TEST_EQ(e * x, x);
    TEST_EQ(n * x, n);
    TEST_EQ(x * y, y * x);
    TEST_EQ((x * y) * z, x * (y * z));

    TEST_EQ(rec(rec(x)), x);
    TEST_EQ(rec(x * y), rec(x) * rec(y));
    TEST_EQ(e / x, rec(x));
    TEST_EQ(x / x, e);

    TEST_EQ(a * (x * y), (a * x) * y);
    TEST_EQ(rec(a * x), rec(a) * rec(x));

    TEST_EQ(abs(x * y), abs(x) * abs(y));
    TEST_EQ(max(e, n), e);
    TEST_EQ(max(-x, -y), -min(x, y));
    TEST_EQ(min(e, n), n);
    TEST_EQ(min(-x, -y), -max(x, y));

    TEST_EQ(sum(n), T(0));
    TEST_EQ(sum(e), T(e.size()));
    TEST_EQ(sum(x + y), sum(x) + sum(y));
    TEST_EQ(sum(a * x), a * sum(x));
    TEST_EQ(maxval(n), T(0));
    TEST_EQ(maxval(e), T(1));
    TEST_EQ(maxval(max(x, y)), max(maxval(x), maxval(y)));
    TEST_EQ(maxval(x + a), maxval(x) + a);
    TEST_EQ(maxval(abs(a) * x), abs(a) * maxval(x));
    TEST_EQ(minval(n), T(0));
    TEST_EQ(minval(e), T(1));
    TEST_EQ(minval(min(x, y)), min(minval(x), minval(y)));
    TEST_EQ(minval(x + a), minval(x) + a);
    TEST_EQ(minval(abs(a) * x), abs(a) * minval(x));
    TEST_EQ(maxval(-a), -minval(a));

    {
      constexpr std::size_t N = n.size();
      std::array<T, 2 * N> arr0;
      for (T &elt : arr0)
        elt = dist(gen);
      std::array<T, 2 * N> arr;
      for (size_t n = 0; n <= N; ++n) {
        // load, store is identity
        arr = arr0;
        const simt<T> x0 = make_memref(arr[n]);
        make_memref(arr[n]) = x0;
        TEST_EQ(arr, arr0);
        // store, load returns stored value
        arr = arr0;
        make_memref(arr[n]) = x;
        const simt<T> x1 = make_memref(arr[n]);
        TEST_EQ(x1, x);
        // double store returns the second stored value
        arr = arr0;
        make_memref(arr[n]) = x;
        make_memref(arr[n]) = y;
        const simt<T> y1 = make_memref(arr[n]);
        TEST_EQ(y1, y);
      }
    }

    {
      constexpr std::size_t N = n.size();
      std::array<T, 2 * N> arr0;
      for (T &elt : arr0)
        elt = dist(gen);
      std::array<T, 2 * N> arr;
      for (size_t n = 0; n <= N; ++n) {
        const auto ptr = make_memptr(&arr[n]);
        // load, store is identity
        arr = arr0;
        const simt<T> x0 = ptr[n];
        ptr[n] = x0;
        TEST_EQ(arr, arr0);
        // store, load returns stored value
        arr = arr0;
        ptr[n] = x;
        const simt<T> x1 = ptr[n];
        TEST_EQ(x1, x);
        // double store returns the second stored value
        arr = arr0;
        ptr[n] = x;
        ptr[n] = y;
        const simt<T> y1 = ptr[n];
        TEST_EQ(y1, y);
      }
    }
  }
}

template <typename T, typename Tptr>
CCTK_ATTRIBUTE_NOINLINE void axpy(const std::ptrdiff_t n, const T a,
                                  const Tptr x, const Tptr y) {
  forall<T>(0, n, [&](const auto i, const auto mask) { y[i] += a * x[i]; });
}

extern "C" void SIMT_Test(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  num_tests = 0;
  num_failed_tests = 0;
  num_inaccurate_tests = 0;

  CCTK_INFO("Testing SIMT::simt<CCTK_REAL4>...");
  run_tests<CCTK_REAL4>();
  CCTK_INFO("Testing SIMT::simt<CCTK_REAL8>...");
  run_tests<CCTK_REAL8>();

  CCTK_VINFO("Overall results:");
  CCTK_VINFO("  tests run:        %5td", num_tests);
  CCTK_VINFO("  tests inaccurate: %5td", num_inaccurate_tests);
  CCTK_VINFO("  tests failed:     %5td", num_failed_tests);

  std::ptrdiff_t n = 1024;
  CCTK_REAL a = 2;
  std::vector<CCTK_REAL> x(n, 3);
  std::vector<CCTK_REAL> y(n, 4);
  axpy(n, a, x.data(), y.data());
  axpy(n, simt<CCTK_REAL>(a), memptr<CCTK_REAL>(x.data()),
       memptr<CCTK_REAL>(y.data()));

  if (num_failed_tests != 0)
    CCTK_ERROR("SIMT self-tests failed");
}

} // namespace SIMT
