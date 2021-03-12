#ifndef IO_SILO_HXX
#define IO_SILO_HXX

#include <fixmath.hxx>
#include <cctk.h>

#ifdef HAVE_CAPABILITY_Silo
namespace CarpetX {

void InputSilo(const cGH *cctkGH);
void OutputSilo(const cGH *cctkGH);

} // namespace CarpetX
#endif

#endif // #ifndef IO_SILO_HXX
