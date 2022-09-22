#include "fillpatch.hxx"

#include <utility>

namespace CarpetX {

using namespace amrex;

namespace {
void do_nothing() {}
std::function<void()> do_nothing2() { return do_nothing; }
} // namespace

std::function<std::function<void()>()> FillPatch_Sync(
    MultiFab &mfab, const Geometry &geom,
    PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t> &bc) {
  mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                           geom.periodicity());
  return [&mfab, &bc]() -> std::function<void()> {
    mfab.FillBoundary_finish();
    bc(mfab, 0, mfab.nComp(), mfab.nGrowVect(), 0.0, 0);
    return do_nothing;
  };
}

std::function<std::function<void()>()> FillPatch_ProlongateGhosts(
    MultiFab &mfab, const MultiFab &cmfab, const Geometry &cgeom,
    const Geometry &fgeom,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &cbc,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &fbc,
    Interpolater *const mapper, const Vector<BCRec> &bcrecs) {
  const IntVect &nghosts = mfab.nGrowVect();
  if (nghosts.max() == 0)
    return do_nothing2;

  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
      mfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);

  // Synchronize
  mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                           fgeom.periodicity());

  if (fpc.ba_crse_patch.empty()) {

    return [&mfab, &fbc]() -> std::function<void()> {
      // Finish synchronizing
      mfab.FillBoundary_finish();
      fbc(mfab, 0, mfab.nComp(), mfab.nGrowVect(), 0.0, 0);
      return do_nothing;
    };

  } else {

    // Copy parts of coarse grid into temporary buffer
    MultiFab *const mfab_crse_patch_ptr =
        new MultiFab(make_mf_crse_patch<MultiFab>(fpc, ncomps));
    MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;
    mf_set_domain_bndry(mfab_crse_patch, cgeom);

    // This is not local
    mfab_crse_patch.ParallelCopy_nowait(cmfab, 0, 0, ncomps, IntVect{0},
                                        mfab_crse_patch.nGrowVect(),
                                        cgeom.periodicity());

    return [&mfab, &cgeom, &fgeom, &cbc, &fbc, mapper, &bcrecs, &fpc,
            mfab_crse_patch_ptr]() -> std::function<void()> {
      const IntVect &nghosts = mfab.nGrowVect();
      const int ncomps = mfab.nComp();
      const IntVect ratio{2, 2, 2};
      MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;

      // Finish synchronizing
      mfab.FillBoundary_finish();
      fbc(mfab, 0, mfab.nComp(), mfab.nGrowVect(), 0.0, 0);

      // Finish copying parts of coarse grid into temporary buffer
      mfab_crse_patch.ParallelCopy_finish();
      cbc(mfab_crse_patch, 0, ncomps, mfab_crse_patch.nGrowVect(), 0.0, 0);

      MultiFab *const mfab_fine_patch_ptr =
          new MultiFab(make_mf_fine_patch<MultiFab>(fpc, ncomps));
      MultiFab &mfab_fine_patch = *mfab_fine_patch_ptr;

      // Interpolate coarse buffer into fine buffer (in space, local)
      FillPatchInterp(mfab_fine_patch, 0, mfab_crse_patch, 0, ncomps,
                      IntVect(0), cgeom, fgeom,
                      grow(convert(fgeom.Domain(), mfab.ixType()), nghosts),
                      ratio, mapper, bcrecs, 0);

      // Copy fine buffer into destination
      mfab.ParallelCopy_nowait(mfab_fine_patch, 0, 0, ncomps, IntVect{0},
                               nghosts);

      delete mfab_crse_patch_ptr;

      return [&mfab, mfab_fine_patch_ptr]() {
        // Finish copying fine buffer into destination
        mfab.ParallelCopy_finish();

        delete mfab_fine_patch_ptr;
      };
    };
  }
}

void FillPatch_Regrid(
    MultiFab &mfab, const MultiFab &cmfab, const MultiFab &fmfab,
    const Geometry &cgeom, const Geometry &fgeom,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &cbc,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &fbc,
    Interpolater *const mapper, const Vector<BCRec> &bcrecs) {
  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const IntVect &nghosts = mfab.nGrowVect();
  const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
      fmfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);

  if (!fpc.ba_crse_patch.empty()) {
    MultiFab mfab_crse_patch = make_mf_crse_patch<MultiFab>(fpc, ncomps);
    mf_set_domain_bndry(mfab_crse_patch, cgeom);

    mfab_crse_patch.ParallelCopy(cmfab, 0, 0, ncomps, IntVect{0}, nghosts,
                                 cgeom.periodicity());
    cbc(mfab_crse_patch, 0, ncomps, nghosts, 0.0, 0);

    MultiFab mfab_fine_patch = make_mf_fine_patch<MultiFab>(fpc, ncomps);

    // In space, local
    FillPatchInterp(mfab_fine_patch, 0, mfab_crse_patch, 0, ncomps, IntVect(0),
                    cgeom, fgeom,
                    grow(convert(fgeom.Domain(), mfab.ixType()), nghosts),
                    ratio, mapper, bcrecs, 0);

    mfab.ParallelCopy_nowait(mfab_fine_patch, 0, 0, ncomps, IntVect{0},
                             nghosts);
    mfab.ParallelCopy_finish();
  }

  mfab.ParallelCopy(fmfab, 0, 0, ncomps, IntVect{0}, nghosts,
                    fgeom.periodicity());
  fbc(mfab, 0, ncomps, nghosts, 0.0, 0);
}

} // namespace CarpetX
