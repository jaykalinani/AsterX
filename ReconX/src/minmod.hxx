#ifndef RECONX_MINMOD_HXX
#define RECONX_MINMOD_HXX

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
inline CCTK_ATTRIBUTE_ALWAYS_INLINE CCTK_DEVICE CCTK_HOST T minmod(const T &x,
                                                                   const T &y) {
  if (signbit(x) != signbit(y))
    return T(0);
  if (fabs(x) < fabs(y))
    return x;
  else
    return y;
}

} // namespace ReconX

#endif // RECONX_MINMOD_HXX
