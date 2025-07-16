#include <loop_device.hxx>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <cmath>

#include "aster_utils.hxx"

namespace AsterAnalysis {
using namespace std;
using namespace Loop;
using namespace EOSX;
using namespace AsterUtils;

extern "C" void AsterAnalysis_MHD_Spherical(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTSX_AsterAnalysis_MHD_Spherical;
  DECLARE_CCTK_PARAMETERS;

}

} // namespace AsterAnalysis
