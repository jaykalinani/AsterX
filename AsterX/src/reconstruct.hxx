#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>

namespace AsterX {

enum class reconstruction_t { Godunov, minmod, monocentral, ppm };

CCTK_DEVICE array<CCTK_REAL, 2>
reconstruct(const GF3D2<const CCTK_REAL> &gf_var, const PointDesc &p,
            reconstruction_t reconstruction, int dir);

} // namespace AsterX
