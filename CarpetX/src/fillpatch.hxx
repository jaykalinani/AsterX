#ifndef FILLPATCH_HXX
#define FILLPATCH_HXX

#include <functional>

#include "driver.hxx"

namespace CarpetX {

// Sync
std::function<std::function<void()>()> FillPatch_Sync(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    amrex::MultiFab &mfab, const amrex::Geometry &geom,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<
        GHExt::PatchData::LevelData::GroupData::apply_physbcs_t> &bc);

// Prolongate (but do not sync) ghosts. Expects coarse mfab synced (but not
// necessarily ghost-prolongated).
std::function<std::function<void()>()> FillPatch_ProlongateGhosts(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &cgeom, const amrex::Geometry &fgeom,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<
        GHExt::PatchData::LevelData::GroupData::apply_physbcs_t> &cbc,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<
        GHExt::PatchData::LevelData::GroupData::apply_physbcs_t> &fbc,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

#warning "TODO: Restrict"

// Prolongate and sync interior. Expects coarse mfab prolongated and synced.
// ("InterpFromCoarseLevel")
void FillPatch_NewLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &cgeom, const amrex::Geometry &fgeom,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<
        GHExt::PatchData::LevelData::GroupData::apply_physbcs_t> &cbc,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<
        GHExt::PatchData::LevelData::GroupData::apply_physbcs_t> &fbc,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

// ("FillPatchTwoLevels")
void FillPatch_RemakeLevel(
    const GHExt::PatchData::LevelData::GroupData &groupdata,
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::MultiFab &fmfab, const amrex::Geometry &cgeom,
    const amrex::Geometry &fgeom,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<
        GHExt::PatchData::LevelData::GroupData::apply_physbcs_t> &cbc,
    GHExt::PatchData::LevelData::GroupData::CactusPhysBCFunct<
        GHExt::PatchData::LevelData::GroupData::apply_physbcs_t> &fbc,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

} // namespace CarpetX

#endif // #ifndef FILLPATCH_HXX
