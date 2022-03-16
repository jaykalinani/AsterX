#include "driver.hxx"
#include "loop_device.hxx"
#include "schedule.hxx"
#include "symmetries.hxx"

#include <vect.hxx>

#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Parameters.h>

#include <array>

namespace CarpetX {

void ApplySymmetries(
    const cGH *restrict const cctkGH,
    const GHExt::PatchData::LevelData::GroupData &restrict groupdata, int tl) {
  DECLARE_CCTK_PARAMETERS;

  // Find faces with reflection symmetry
  const std::array<std::array<bool, 3>, 2> is_reflect{{
      {{bool(reflection_x), bool(reflection_y), bool(reflection_z)}},
      {{bool(reflection_upper_x), bool(reflection_upper_y),
        bool(reflection_upper_z)}},
  }};
  // If there are no faces then we are already done
  bool any_reflection = false;
  for (int d = 0; d < 3; ++d)
    for (int f = 0; f < 2; ++f)
      any_reflection |= is_reflect[f][d];
  if (!any_reflection)
    return;

  // Find which of the 27 sub-boxes need to be filled by the symmetry
  std::array<std::array<std::array<bool, 3>, 3>, 3> there;
  for (int k = -1; k <= +1; ++k)
    for (int j = -1; j <= +1; ++j)
      for (int i = -1; i <= +1; ++i)
        there[i + 1][j + 1][k + 1] = false;
  for (int d = 0; d < 3; ++d) {
    for (int f = 0; f < 2; ++f) {
      if (is_reflect[f][d]) {
        // Select this face
        for (int k = -1; k <= +1; ++k) {
          for (int j = -1; j <= +1; ++j) {
            for (int i = -1; i <= +1; ++i) {
              const std::array<int, 3> idx{i, j, k};
              if (idx[d] == (f == 0 ? -1 : +1))
                there[i + 1][j + 1][k + 1] = true;
            }
          }
        }
        //
      }
    }
  }

  // Loop over all blocks on the current level of the current patch
  loop_over_blocks(cctkGH,
                   active_levels_t(groupdata.level, groupdata.level + 1,
                                   groupdata.patch, groupdata.patch + 1),
                   [&](const int patch, const int level, const int index,
                       const int block, const cGH *restrict const cctkGH) {
                     // Get memory layout for this grid function group
                     const Loop::GF3D2layout layout(cctkGH, groupdata.indextype,
                                                    groupdata.nghostzones);
                     const Loop::GridDescBaseDevice grid(cctkGH);

                     for (int vi = 0; vi < groupdata.numvars; ++vi) {
                       CCTK_REAL *const ptr =
                           static_cast<CCTK_REAL *>(CCTK_VarDataPtrI(
                               cctkGH, tl, groupdata.firstvarindex + vi));
                       const Loop::GF3D2<CCTK_REAL> var(layout, ptr);
                       const Arith::vect<int, dim> CI = groupdata.indextype;
                       grid.loop_there_device_idx(
                           groupdata.indextype, there, groupdata.nghostzones,
                           [=] CCTK_DEVICE(const Loop::PointDesc &p)
                               CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                 const auto DI = p.I - p.I0;
                                 const auto J = p.I0 - DI - (1 - CI) * p.NI;
                                 var(p.I) = var(J);
                               });
                     } // for vi
                   });
}

} // namespace CarpetX
