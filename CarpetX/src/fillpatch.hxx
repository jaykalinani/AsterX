#ifndef FILLPATCH_HXX
#define FILLPATCH_HXX

#include <functional>

#include "driver.hxx"

namespace CarpetX {

// Sync
std::function<std::function<void()>()> FillPatch_Sync(
    amrex::MultiFab &mfab, const amrex::Geometry &geom,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &bc);

// Prolongate (but do not sync) ghosts. Expects coarse mfab synced (but not
// necessarily ghost-prolongated).
std::function<std::function<void()>()> FillPatch_ProlongateGhosts(
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::Geometry &cgeom, const amrex::Geometry &fgeom,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &cbc,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &fbc,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

#warning "TODO: Restrict"

// Prolongate and sync interior. Expects coarse mfab prolongated and synced.
void FillPatch_Regrid(
    amrex::MultiFab &mfab, const amrex::MultiFab &cmfab,
    const amrex::MultiFab &fmfab, const amrex::Geometry &cgeom,
    const amrex::Geometry &fgeom,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &cbc,
    amrex::PhysBCFunct<GHExt::PatchData::LevelData::GroupData::apply_physbcs_t>
        &fbc,
    amrex::Interpolater *mapper, const amrex::Vector<amrex::BCRec> &bcrecs);

} // namespace CarpetX

#endif // #ifndef FILLPATCH_HXX
