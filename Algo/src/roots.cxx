#include "roots.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>

#include <cassert>
#include <cmath>
#include <limits>
#include <tuple>

namespace Algo {

extern "C" void Test_roots(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;

  auto fn = [](auto x) { return x * x - 2; };
  auto fnd = [](auto x) { return std::make_tuple(x * x - 2, 2 * x); };

  {
    const int minbits = std::numeric_limits<double>::digits - 4;
    const int maxiters = 100;
    int iters;
    auto [lo, hi] = bisect(fn, 1.0, 2.0, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    assert(hi >= lo && hi - lo <= std::scalbn(2.0, -minbits));
    assert(fn(lo) <= 0 && fn(hi) >= 0);
    CCTK_VINFO("Test_bisect succeeded in %d iterations", iters);
  }

  {
    const int minbits = std::numeric_limits<double>::digits - 4;
    const int maxiters = 100;
    int iters;
    auto [lo, hi] =
        bracket_and_solve_root(fn, 1.0, 2.0, true, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    assert(hi >= lo && hi - lo <= std::scalbn(2.0, -minbits));
    assert(fn(lo) <= 0 && fn(hi) >= 0);
    CCTK_VINFO("Test_bracket_and_solve_root succeeded in %d iterations", iters);
  }

  {
    const int minbits = std::numeric_limits<double>::digits - 4;
    const int maxiters = 100;
    int iters;
    auto [lo, hi] = brent(fn, 1.0, 2.0, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    assert(hi >= lo && hi - lo <= std::scalbn(2.0, -minbits));
    assert(fn(lo) <= 0 && fn(hi) >= 0);
    CCTK_VINFO("Test_brent succeeded in %d iterations", iters);
  }

  {
    const int minbits =
        static_cast<int>(0.6 * std::numeric_limits<double>::digits);
    const int maxiters = 100;
    int iters;
    auto x = newton_raphson(fnd, 1.0, 0.0, 10.0, minbits, maxiters, iters);
    // CCTK_VINFO("maxiters=%d iters=%d", maxiters, iters);
    // CCTK_VINFO("lo=%.17g hi=%.17g", double(lo), double(hi));
    assert(iters < maxiters);
    double delta = std::scalbn(1.0, -minbits);
    assert(fn(x - delta) * fn(x + delta) < 0);
    CCTK_VINFO("Test_newton_raphson succeeded in %d iterations", iters);
  }
}

} // namespace Algo
