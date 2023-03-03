#ifndef IO_HXX
#define IO_HXX

#include <cctk.h>

namespace CarpetX {

void RecoverGridStructure(cGH *cctkGH);
void RecoverGH(const cGH *cctkGH);
void InputGH(const cGH *cctkGH);

int OutputGH(const cGH *cctkGH);

} // namespace CarpetX

#endif // #ifndef IO_HXX
