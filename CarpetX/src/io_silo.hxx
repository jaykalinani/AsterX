#ifndef IO_SILO_HXX
#define IO_SILO_HXX

#include <cctk.h>

#ifdef HAVE_CAPABILITY_Silo
namespace CarpetX {

void OutputSilo(const cGH *restrict const cctkGH);

}
#endif

#endif // #ifndef IO_SILO_HXX
