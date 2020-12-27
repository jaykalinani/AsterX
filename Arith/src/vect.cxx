#include "vect.hxx"

namespace Arith {
using namespace std;

namespace TestVect {
constexpr vect<int, 3> vect1{0, 1, 2};
static_assert(vect1[0] == 0, "");
static_assert(vect1[1] == 1, "");
static_assert(vect1[2] == 2, "");

constexpr vect<vect<int, 3>, 3> vect2{
    {100, 101, 102},
    {110, 111, 112},
    {120, 121, 122},
};
static_assert(vect2[0][0] == 100, "");
static_assert(vect2[0][1] == 101, "");
static_assert(vect2[0][2] == 102, "");
static_assert(vect2[1][0] == 110, "");
static_assert(vect2[1][1] == 111, "");
static_assert(vect2[1][2] == 112, "");
static_assert(vect2[2][0] == 120, "");
static_assert(vect2[2][1] == 121, "");
static_assert(vect2[2][2] == 122, "");

} // namespace TestVect

} // namespace Arith
