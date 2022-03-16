#ifndef SYMMETRIES_HXX
#define SYMMETRIES_HXX

#include "driver.hxx"

#include <fixmath.hxx>
#include <cctk.h>

namespace CarpetX {

void ApplySymmetries(
    const cGH *restrict const cctkGH,
    const GHExt::PatchData::LevelData::GroupData &restrict groupdata, int tl);

}

#endif // #ifndef SYMMETRIES_HXX
