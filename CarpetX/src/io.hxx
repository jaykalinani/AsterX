#ifndef IO_HXX
#define IO_HXX

#include <cctk.h>

namespace CarpetX {

int InputGH(const cGH *cctkGH);
int OutputGH(const cGH *cctkGH);

} // namespace CarpetX

#endif // #ifndef IO_HXX
