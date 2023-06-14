#ifndef RECONX_MONOCENTRAL_HXX
#define RECONX_MONOCENTRAL_HXX

#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <cassert>

#include "reconx_utils.hxx"

namespace ReconX {

using namespace std;
using namespace Arith;

template <typename T>
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T
monocentral(const T &x, const T &y) {
  if (sgn(x) != sgn(y))
    return 0;
  else
    return sgn(x) * min(2 * fabs(x), min(2 * fabs(y), fabs(x + y) / 2));
}

} // namespace ReconX

#endif // RECONX_MINMOD_HXX
