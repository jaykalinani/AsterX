#ifndef DEFS_HH
#define DEFS_HH

#include <simd.hxx>

#include <array>
#include <cmath>
#include <tuple>

namespace Hydro {
using namespace Arith;
using namespace std;

constexpr int dim = 3;

typedef simd<CCTK_REAL> CCTK_REALVEC;
typedef simdl<CCTK_REAL> CCTK_BOOLVEC;
constexpr int vsize = tuple_size_v<CCTK_REALVEC>;

} // namespace Hydro

#endif // #ifndef DEFS_HH
