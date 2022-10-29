#include "fillpatch.hxx"

#include <utility>

#include <AMReX_FillPatchUtil.H>

namespace CarpetX {

using namespace amrex;

using apply_physbcs_t = GHExt::PatchData::LevelData::GroupData::apply_physbcs_t;

namespace {
void do_nothing() {}
std::function<void()> do_nothing2() { return do_nothing; }
} // namespace

std::function<std::function<void()>()> FillPatch_Sync(
    const GHExt::PatchData::LevelData::GroupData &groupdata, MultiFab &mfab,
    const Geometry &geom,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<apply_physbcs_t>
        &bc) {
  mfab.FillBoundary_nowait(0, mfab.nComp(), mfab.nGrowVect(),
                           geom.periodicity());
  return [&groupdata, &mfab, &bc]() -> std::function<void()> {
    mfab.FillBoundary_finish();
    bc(groupdata, mfab, 0, mfab.nComp(), mfab.nGrowVect(), 0.0, 0);
    return do_nothing;
  };
}

std::function<std::function<void()>()> FillPatch_ProlongateGhosts(
    const GHExt::PatchData::LevelData::GroupData &groupdata, MultiFab &mfab,
    const MultiFab &cmfab, const Geometry &cgeom, const Geometry &fgeom,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<apply_physbcs_t>
        &cbc,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<apply_physbcs_t>
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

    return [&groupdata, &mfab, &fbc]() -> std::function<void()> {
      // Finish synchronizing
      mfab.FillBoundary_finish();
      fbc(groupdata, mfab, 0, mfab.nComp(), mfab.nGrowVect(), 0.0, 0);
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

    return [&groupdata, &mfab, &cgeom, &fgeom, &cbc, &fbc, mapper, &bcrecs,
            &fpc, mfab_crse_patch_ptr]() -> std::function<void()> {
      const IntVect &nghosts = mfab.nGrowVect();
      const int ncomps = mfab.nComp();
      const IntVect ratio{2, 2, 2};
      MultiFab &mfab_crse_patch = *mfab_crse_patch_ptr;

      // Finish synchronizing
      mfab.FillBoundary_finish();
      fbc(groupdata, mfab, 0, mfab.nComp(), mfab.nGrowVect(), 0.0, 0);

      // Finish copying parts of coarse grid into temporary buffer
      mfab_crse_patch.ParallelCopy_finish();
      cbc(groupdata, mfab_crse_patch, 0, ncomps, mfab_crse_patch.nGrowVect(),
          0.0, 0);

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

void FillPatch_NewLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata, MultiFab &mfab,
    const MultiFab &cmfab, const Geometry &cgeom, const Geometry &fgeom,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<apply_physbcs_t>
        &cbc,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<apply_physbcs_t>
        &fbc,
    Interpolater *const mapper, const Vector<BCRec> &bcrecs) {

  const int ncomps = mfab.nComp();
  const IntVect ratio{2, 2, 2};
  const IntVect &nghosts = mfab.nGrowVect();
  // const EB2::IndexSpace *const index_space = nullptr;

  const InterpolaterBoxCoarsener &coarsener = mapper->BoxCoarsener(ratio);

  const BoxArray &ba = mfab.boxArray();
  const DistributionMapping &dm = mfab.DistributionMap();

  const IndexType &ixtype = ba.ixType();
  assert(ixtype == cmfab.boxArray().ixType());

  // const FabArrayBase::FPinfo &fpc = FabArrayBase::TheFPinfo(
  //     fmfab, mfab, nghosts, coarsener, fgeom, cgeom, index_space);
  //
  // if (!fpc.ba_crse_patch.empty()) {
  //   MultiFab mfab_crse_patch = make_mf_crse_patch<MultiFab>(fpc, ncomps);
  //   mf_set_domain_bndry(mfab_crse_patch, cgeom);
  //
  //   mfab_crse_patch.ParallelCopy(cmfab, 0, 0, ncomps, IntVect{0}, nghosts,
  //                                cgeom.periodicity());
  //   cbc(groupdata, mfab_crse_patch, 0, ncomps, nghosts, 0.0, 0);
  //
  //   MultiFab mfab_fine_patch = make_mf_fine_patch<MultiFab>(fpc, ncomps);
  //
  //   // In space, local
  //   FillPatchInterp(mfab_fine_patch, 0, mfab_crse_patch, 0, ncomps,
  //   IntVect(0),
  //                   cgeom, fgeom,
  //                   grow(convert(fgeom.Domain(), mfab.ixType()), nghosts),
  //                   ratio, mapper, bcrecs, 0);
  //
  //   mfab.ParallelCopy_nowait(mfab_fine_patch, 0, 0, ncomps, IntVect{0},
  //                            nghosts);
  //   mfab.ParallelCopy_finish();
  // }
  //
  // mfab.ParallelCopy(fmfab, 0, 0, ncomps, IntVect{0}, nghosts,
  //                   fgeom.periodicity());
  // fbc(groupdata, mfab, 0, ncomps, nghosts, 0.0, 0);

  // vvv new

  // Suffix `_g` is for "with ghosts added"
  Box fdomain_g(amrex::convert(fgeom.Domain(), mfab.ixType()));
  for (int d = 0; d < dim; ++d)
    if (fgeom.isPeriodic(d))
      fdomain_g.grow(d, nghosts[d]);

  const int nboxes = ba.size();
  BoxArray cba_g(nboxes);
  for (int i = 0; i < nboxes; ++i) {
    Box box = amrex::convert(amrex::grow(ba[i], nghosts), ixtype);
    box &= fdomain_g;
    cba_g.set(i, coarsener.doit(box));
  }
  MultiFab cmfab_g(cba_g, dm, ncomps, 0);
  mf_set_domain_bndry(cmfab_g, cgeom);

  cmfab_g.ParallelCopy(cmfab, 0, 0, ncomps, cgeom.periodicity());

  cbc(groupdata, cmfab_g, 0, ncomps, cmfab_g.nGrowVect(), 0.0, 0);

  FillPatchInterp(mfab, 0, cmfab_g, 0, ncomps, nghosts, cgeom, fgeom, fdomain_g,
                  ratio, mapper, bcrecs, 0);

  fbc(groupdata, mfab, 0, ncomps, nghosts, 0.0, 0);
}

void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata, MultiFab &mfab,
    const MultiFab &cmfab, const MultiFab &fmfab, const Geometry &cgeom,
    const Geometry &fgeom,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<apply_physbcs_t>
        &cbc,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<apply_physbcs_t>
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
    cbc(groupdata, mfab_crse_patch, 0, ncomps, nghosts, 0.0, 0);

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
  fbc(groupdata, mfab, 0, ncomps, nghosts, 0.0, 0);
}

} // namespace CarpetX
