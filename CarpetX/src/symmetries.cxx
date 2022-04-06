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
    const GHExt::PatchData::LevelData::GroupData &restrict groupdata,
    const int tl) {
  DECLARE_CCTK_PARAMETERS;

  const auto &symmetries = groupdata.leveldata().patchdata().symmetries;

  // Find which of the 27 sub-boxes need to be filled by which
  // symmetry. Periodicity is already handled by AMReX. If there are
  // multiple symmetries applicable, choose one of them.
  std::array<std::array<std::array<bool, 3>, 3>, 3> there_dirichlet,
      there_reflection;
  for (int k = -1; k <= +1; ++k) {
    for (int j = -1; j <= +1; ++j) {
      for (int i = -1; i <= +1; ++i) {
        there_dirichlet[i + 1][j + 1][k + 1] = false;
        there_reflection[i + 1][j + 1][k + 1] = false;
      }
    }
  }
  for (int d = 0; d < 3; ++d) {
    for (int f = 0; f < 2; ++f) {
      // Select this face
      for (int k = -1; k <= +1; ++k) {
        for (int j = -1; j <= +1; ++j) {
          for (int i = -1; i <= +1; ++i) {
            const std::array<int, 3> idx{i, j, k};
            if (idx[d] == (f == 0 ? -1 : +1)) {
              switch (symmetries[f][d]) {
              case symmetry_t::dirichlet:
                there_dirichlet[i + 1][j + 1][k + 1] = true;
                break;
              case symmetry_t::reflection:
                there_reflection[i + 1][j + 1][k + 1] = true;
                break;
              default:
                // do nothing
                break;
              }
            }
          }
        }
      }
      //
    }
  }

  for (int k = -1; k <= +1; ++k)
    for (int j = -1; j <= +1; ++j)
      for (int i = -1; i <= +1; ++i)
        there_reflection[i + 1][j + 1][k + 1] &=
            !there_dirichlet[i + 1][j + 1][k + 1];
  for (int k = -1; k <= +1; ++k)
    for (int j = -1; j <= +1; ++j)
      for (int i = -1; i <= +1; ++i)
        assert(there_reflection[i + 1][j + 1][k + 1] +
                   there_dirichlet[i + 1][j + 1][k + 1] <=
               1);
  bool any_dirichlet = false, any_reflection = false;
  for (int k = -1; k <= +1; ++k) {
    for (int j = -1; j <= +1; ++j) {
      for (int i = -1; i <= +1; ++i) {
        any_dirichlet |= there_dirichlet[i + 1][j + 1][k + 1];
        any_reflection |= there_reflection[i + 1][j + 1][k + 1];
      }
    }
  }

  // If there are no faces then we are already done
  if (!any_dirichlet && !any_reflection)
    return;

#warning                                                                       \
    "TODO: check that inputs are valid. set outputs as valid afterwards (?) shouldn't that be handled by the caller, e.g. sync?"

  // Loop over all blocks on the current level of the current patch
  loop_over_blocks(
      active_levels_t(groupdata.level, groupdata.level + 1, groupdata.patch,
                      groupdata.patch + 1),
      [&](const int patch, const int level, const int index, const int block,
          const cGH *restrict const cctkGH) {
        // Get memory layout for this grid function group
        const Loop::GF3D2layout layout(cctkGH, groupdata.indextype,
                                       groupdata.nghostzones);
        const Loop::GridDescBaseDevice grid(cctkGH);
        const Arith::vect<int, dim> CI = groupdata.indextype;

        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          CCTK_REAL *const ptr = static_cast<CCTK_REAL *>(
              CCTK_VarDataPtrI(cctkGH, tl, groupdata.firstvarindex + vi));
          const Loop::GF3D2<CCTK_REAL> var(layout, ptr);
          const CCTK_REAL dirichlet_value = groupdata.dirichlet_values.at(vi);
          const Arith::vect<int, dim> P = groupdata.parities.at(vi);

#warning "TODO: mask with cctkGH->cctk_bbox"

          if (any_dirichlet)
            grid.loop_there_device_idx(groupdata.indextype, there_dirichlet,
                                       groupdata.nghostzones,
                                       [=] CCTK_DEVICE(const Loop::PointDesc &p)
                                           CCTK_ATTRIBUTE_ALWAYS_INLINE {
                                             var(p.I) = dirichlet_value;
                                           });

          if (any_reflection)
            grid.loop_there_device_idx(
                groupdata.indextype, there_reflection, groupdata.nghostzones,
                [=] CCTK_DEVICE(const Loop::PointDesc &p)
                    CCTK_ATTRIBUTE_ALWAYS_INLINE {
                      int parity = 1;
                      for (int d = 0; d < dim; ++d)
                        if (p.NI[d] != 0)
                          parity *= P[d];
                      const auto DI = p.I0 - p.I;
                      const auto J = p.I0 + DI + CI * p.NI;
                      var(p.I) = parity * var(J);
                    });

        } // for vi
      });
}

} // namespace CarpetX
