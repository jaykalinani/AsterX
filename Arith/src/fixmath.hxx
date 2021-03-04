#ifndef FIXMATH_HXX
#define FIXMATH_HXX

// Include this file BEFORE including <cctk.h>, best before including
// any other include files

// This provides broken `#defines` for `isnan` etc.
#include <cctk.h>

// Disable the broken definitions
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#endif // #ifndef FIXMATH_HXX
