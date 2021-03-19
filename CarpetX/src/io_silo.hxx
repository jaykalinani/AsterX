#ifndef IO_SILO_HXX
#define IO_SILO_HXX

#include <fixmath.hxx>
#include <cctk.h>

#ifdef HAVE_CAPABILITY_Silo
namespace CarpetX {

void InputSilo(const cGH *cctkGH);
enum class output_type_t { scheduled, checkpoint };
void OutputSilo(const cGH *cctkGH, output_type_t output_type);

} // namespace CarpetX
#endif

#endif // #ifndef IO_SILO_HXX
