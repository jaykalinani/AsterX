#ifndef IO_HXX
#define IO_HXX

#include <cctk.h>
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

namespace CarpetX {

int OutputGH(const cGH *cctkGH);
} // namespace CarpetX

#endif // #ifndef IO_HXX
