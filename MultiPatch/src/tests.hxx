#ifndef MULTIPATCH_TESTS_HXX
#define MULTIPATCH_TESTS_HXX

#include <cmath>
#include <string>
#include <random>

namespace MultiPatchTests {

constexpr const int random_seed = 100;

/**
 * Test two floating point values for approximate equality.
 *
 * This function is based on
 * https://docs.julialang.org/en/v1/base/math/#Base.isapprox
 *
 * @param x The first number to compare.
 * @param y The second number to compare.
 * @param atol The absolute tolerance for the comparison.
 * @param rtol The relative tolerance of the comparison.
 * @return True if x ~ y, false otherwise.
 */
template <typename fp_type>
bool isapprox(fp_type x, fp_type y, fp_type atol = 0.0) {
  using std::abs;
  using std::max;
  using std::sqrt;

  const fp_type rtol =
      atol > 0.0 ? 0.0 : sqrt(std::numeric_limits<fp_type>::epsilon());
  return abs(x - y) <= max(atol, rtol * max(abs(x), abs(y)));
}

enum class string_color { none, green, red };

template <string_color color>
constexpr const std::string colored(const std ::string &string) {
  if constexpr (color == string_color::red) {
    std::string msg{"\033[31;1m"};
    msg += string;
    msg += "\033[0m";
    return msg;

  } else if constexpr (color == string_color::green) {
    std::string msg{"\033[32;1m"};
    msg += string;
    msg += "\033[0m";

    return msg;
  }

  return string;
}

} // namespace MultiPatchTests

/**
 * Runs all Cartesian patch tests.
 *
 * This function needs to be scheduled scheduled at the wragh bin
 */
extern "C" void run_cartesian_tests();

/**
 * Runs all Cake patch tests.
 *
 * This function needs to be scheduled scheduled at the wragh bin
 */
extern "C" void run_cake_tests();

#endif // MULTIPATCH_TESTS_HXX