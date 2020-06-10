#ifndef DRIVER_HXX
#define DRIVER_HXX

#include "loop.hxx"
#include "valid.hxx"

#include <rational.hxx>

#include <cctk.h>

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace CarpetX {
using namespace std;
using namespace Arith;

using Loop::dim;

using rat64 = rational<int64_t>;

// Taken from
// <https://stackoverflow.com/questions/27440953/stdunique-ptr-for-c-functions-that-need-free>
struct free_deleter {
  template <typename T> void operator()(T *p) const {
    std::free(const_cast<std::remove_const_t<T> *>(p));
  }
};
template <typename T> using unique_C_ptr = std::unique_ptr<T, free_deleter>;
static_assert(sizeof(char *) == sizeof(unique_C_ptr<char>),
              ""); // ensure no overhead

static_assert(AMREX_SPACEDIM == dim,
              "AMReX's AMREX_SPACEDIM must be the same as Cactus's cctk_dim");

static_assert(is_same<amrex::Real, CCTK_REAL>::value,
              "AMReX's Real type must be the same as Cactus's CCTK_REAL");

class CactusAmrCore final : public amrex::AmrCore {
public:
  CactusAmrCore();
  CactusAmrCore(const amrex::RealBox *rb, int max_level_in,
                const amrex::Vector<int> &n_cell_in, int coord = -1,
                amrex::Vector<amrex::IntVect> ref_ratios =
                    amrex::Vector<amrex::IntVect>(),
                const int *is_per = nullptr);
  CactusAmrCore(const amrex::RealBox &rb, int max_level_in,
                const amrex::Vector<int> &n_cell_in, int coord,
                amrex::Vector<amrex::IntVect> const &ref_ratios,
                amrex::Array<int, AMREX_SPACEDIM> const &is_per);
  CactusAmrCore(const amrex::AmrCore &rhs) = delete;
  CactusAmrCore &operator=(const amrex::AmrCore &rhs) = delete;

  virtual ~CactusAmrCore() override;

  virtual void ErrorEst(int level, amrex::TagBoxArray &tags, amrex::Real time,
                        int ngrow) override;
  virtual void
  MakeNewLevelFromScratch(int level, amrex::Real time,
                          const amrex::BoxArray &ba,
                          const amrex::DistributionMapping &dm) override;
  virtual void
  MakeNewLevelFromCoarse(int level, amrex::Real time, const amrex::BoxArray &ba,
                         const amrex::DistributionMapping &dm) override;
  virtual void RemakeLevel(int level, amrex::Real time,
                           const amrex::BoxArray &ba,
                           const amrex::DistributionMapping &dm) override;
  virtual void ClearLevel(int level) override;
};

// Cactus grid hierarchy extension
struct GHExt {

  // AMReX grid structure
  // TODO: Remove unique_ptr once amrex::AmrCore has move constructors
  unique_ptr<CactusAmrCore> amrcore;

  struct CommonGroupData {
    int groupindex;
    int firstvarindex;
    int numvars;

    bool do_checkpoint; // whether to checkpoint
    bool do_restrict;   // whether to restrict

    vector<vector<why_valid_t> > valid; // [time level][var index]

    // TODO: add poison_invalid and check_valid functions
  };

  struct GlobalData {
    // all data that exists on all levels

    // For subcycling in time, there really should be one copy of each
    // integrated grid scalar per level. We don't do that yet; insted,
    // we assume that grid scalars only hold "analysis" data.

    struct ScalarGroupData : public CommonGroupData {
      vector<vector<CCTK_REAL> > data; // [time level][var index]
    };
    // TODO: right now this is sized for the total number of groups
    vector<unique_ptr<ScalarGroupData> > scalargroupdata; // [group index]
  };
  GlobalData globaldata;

  struct LevelData {
    int level;

    // This level uses subcycling with respect to the next coarser
    // level. (Ignored for the coarsest level.)
    bool is_subcycling_level;

    // Iteration and time at which this cycle level is valid
    rat64 iteration, delta_iteration;

    // Fabamrex::ArrayBase object holding a cell-centred BoxArray for
    // iterating over grid functions. This stores the grid structure
    // and its distribution over all processes, but holds no data.
    unique_ptr<amrex::FabArrayBase> fab;

    struct GroupData : public CommonGroupData {
      int level;
      const LevelData &leveldata() const;
      LevelData &leveldata();

      array<int, dim> indextype;
      array<int, dim> nghostzones;
      vector<array<int, dim> > parities;

      // each amrex::MultiFab has numvars components
      vector<unique_ptr<amrex::MultiFab> > mfab; // [time level]

      // flux register between this and the next coarser level
      unique_ptr<amrex::FluxRegister> freg;
      // associated flux group indices
      array<int, dim> fluxes; // [dir]
    };
    // TODO: right now this is sized for the total number of groups
    vector<unique_ptr<GroupData> > groupdata; // [group index]
  };
  vector<LevelData> leveldata; // [reflevel]
};

extern unique_ptr<GHExt> ghext;

inline const GHExt::LevelData &GHExt::LevelData::GroupData::leveldata() const {
  return ghext->leveldata.at(level);
}
inline GHExt::LevelData &GHExt::LevelData::GroupData::leveldata() {
  return ghext->leveldata.at(level);
}

amrex::Interpolater *get_interpolator(const array<int, dim> indextype);

typedef void apply_physbcs_t(const amrex::Box &, const amrex::FArrayBox &, int,
                             int, const amrex::Geometry &, CCTK_REAL,
                             const amrex::Vector<amrex::BCRec> &, int, int);
typedef amrex::PhysBCFunct<apply_physbcs_t *> CarpetXPhysBCFunct;
tuple<CarpetXPhysBCFunct, amrex::Vector<amrex::BCRec> >
get_boundaries(const GHExt::LevelData::GroupData &groupdata);

} // namespace CarpetX

#endif // #ifndef DRIVER_HXX
