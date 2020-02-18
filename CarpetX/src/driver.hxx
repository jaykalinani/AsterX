#ifndef DRIVER_HXX
#define DRIVER_HXX

#include "vect.hxx"

#include <cctk.h>
#undef copysign
#undef fpclassify
#undef isfinite
#undef isinf
#undef isnan
#undef isnormal
#undef signbit

#include <AMReX.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Interpolater.H>
#include <AMReX_MultiFab.H>

#include <algorithm>
#include <array>
#include <functional>
#include <memory>
#include <ostream>
#include <tuple>
#include <type_traits>
#include <vector>

namespace CarpetX {
using namespace amrex;
using namespace std;

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

static_assert(is_same<Real, CCTK_REAL>::value,
              "AMReX's Real type must be the same as Cactus's CCTK_REAL");

class CactusAmrCore final : public AmrCore {
public:
  CactusAmrCore();
  CactusAmrCore(const RealBox *rb, int max_level_in,
                const Vector<int> &n_cell_in, int coord = -1,
                Vector<IntVect> ref_ratios = Vector<IntVect>(),
                const int *is_per = nullptr);
  CactusAmrCore(const RealBox &rb, int max_level_in,
                const Vector<int> &n_cell_in, int coord,
                Vector<IntVect> const &ref_ratios,
                Array<int, AMREX_SPACEDIM> const &is_per);
  CactusAmrCore(const AmrCore &rhs) = delete;
  CactusAmrCore &operator=(const AmrCore &rhs) = delete;

  virtual ~CactusAmrCore() override;

  virtual void ErrorEst(int level, TagBoxArray &tags, Real time,
                        int ngrow) override;
  virtual void MakeNewLevelFromScratch(int level, Real time, const BoxArray &ba,
                                       const DistributionMapping &dm) override;
  virtual void MakeNewLevelFromCoarse(int level, Real time, const BoxArray &ba,
                                      const DistributionMapping &dm) override;
  virtual void RemakeLevel(int level, Real time, const BoxArray &ba,
                           const DistributionMapping &dm) override;
  virtual void ClearLevel(int level) override;
};

struct valid_t {
  bool valid_int, valid_outer, valid_ghosts;
  constexpr valid_t() : valid_t(false) {}

  explicit constexpr valid_t(bool b)
      : valid_int(b), valid_outer(b), valid_ghosts(b) {}
  friend constexpr valid_t operator~(const valid_t &x) {
    valid_t r;
    r.valid_int = !x.valid_int;
    r.valid_outer = !x.valid_outer;
    r.valid_ghosts = !x.valid_ghosts;
    return r;
  }
  friend constexpr valid_t operator&(const valid_t &x, const valid_t &y) {
    valid_t r;
    r.valid_int = x.valid_int && y.valid_int;
    r.valid_outer = x.valid_outer && y.valid_outer;
    r.valid_ghosts = x.valid_ghosts && y.valid_ghosts;
    return r;
  }
  friend constexpr valid_t operator|(const valid_t &x, const valid_t &y) {
    valid_t r;
    r.valid_int = x.valid_int || y.valid_int;
    r.valid_outer = x.valid_outer || y.valid_outer;
    r.valid_ghosts = x.valid_ghosts || y.valid_ghosts;
    return r;
  }
  valid_t &operator&=(const valid_t &x) { return *this = *this & x; }
  valid_t &operator|=(const valid_t &x) { return *this = *this | x; }

  constexpr bool valid_all() const {
    return valid_int && valid_outer && valid_ghosts;
  }
  constexpr bool valid_any() const {
    return valid_int || valid_outer || valid_ghosts;
  }

  friend bool operator==(const valid_t &x, const valid_t &y) {
    return make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) ==
           make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
  friend bool operator<(const valid_t &x, const valid_t &y) {
    return make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) <
           make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }

  friend ostream &operator<<(ostream &os, const valid_t v) {
    auto str = [](bool v) { return v ? "VAL" : "INV"; };
    return os << "[int:" << str(v.valid_int) << ",outer:" << str(v.valid_outer)
              << ",ghosts:" << str(v.valid_ghosts) << "]";
  }
  operator string() const {
    ostringstream buf;
    buf << *this;
    return buf.str();
  }
};

inline constexpr valid_t make_valid_int() {
  valid_t valid;
  valid.valid_int = true;
  return valid;
}
inline constexpr valid_t make_valid_outer() {
  valid_t valid;
  valid.valid_outer = true;
  return valid;
}
inline constexpr valid_t make_valid_ghosts() {
  valid_t valid;
  valid.valid_ghosts = true;
  return valid;
}
inline constexpr valid_t make_valid_all() { return ~valid_t(); }

} // namespace CarpetX
namespace std {
using namespace CarpetX;
template <> struct equal_to<valid_t> {
  constexpr bool operator()(const valid_t &x, const valid_t &y) const {
    return make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) ==
           make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
};
template <> struct less<valid_t> {
  constexpr bool operator()(const valid_t &x, const valid_t &y) const {
    return make_tuple(x.valid_int, x.valid_outer, x.valid_ghosts) <
           make_tuple(y.valid_int, y.valid_outer, y.valid_ghosts);
  }
};
} // namespace std
namespace CarpetX {

class why_valid_t {
  valid_t valid;
  function<string()> why_int, why_outer, why_ghosts;

public:
  // The constructor that doesn't give a reason should never be called
  why_valid_t() = delete;
  why_valid_t(const function<string()> &why) : why_valid_t(false, why) {}
  why_valid_t(bool b, const function<string()> &why)
      : why_valid_t(valid_t(b), why) {}
  why_valid_t(const valid_t &val, const function<string()> &why)
      : valid(val), why_int(why), why_outer(why), why_ghosts(why) {}

  const valid_t &get() const { return valid; }

  void set(const valid_t &val, const function<string()> &why) {
    auto old_valid = valid;
    valid = val;
    auto new_valid = valid;
    if (new_valid.valid_int != old_valid.valid_int)
      why_int = why;
    if (new_valid.valid_outer != old_valid.valid_outer)
      why_outer = why;
    if (new_valid.valid_ghosts != old_valid.valid_ghosts)
      why_ghosts = why;
  }
  void set_int(bool b, const function<string()> &why) {
    auto val = valid;
    val.valid_int = b;
    set(val, why);
  }
  void set_outer(bool b, const function<string()> &why) {
    auto val = valid;
    val.valid_outer = b;
    set(val, why);
  }
  void set_ghosts(bool b, const function<string()> &why) {
    auto val = valid;
    val.valid_ghosts = b;
    set(val, why);
  }
  void set_and(const valid_t &val, const function<string()> &why) {
    set(valid & val, why);
  }
  void set_or(const valid_t &val, const function<string()> &why) {
    set(valid | val, why);
  }

  friend ostream &operator<<(ostream &os, const why_valid_t why) {
    return os << why.valid << ","
              << "why{int:" << why.why_int() << ","
              << "outer:" << why.why_outer() << ","
              << "ghosts:" << why.why_ghosts() << "}";
  }
  operator string() const {
    ostringstream buf;
    buf << *this;
    return buf.str();
  }
};

// Cactus grid hierarchy extension
struct GHExt {

  // AMReX grid structure
  // TODO: Remove unique_ptr once AmrCore has move constructors
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

    struct ScalarGroupData : public CommonGroupData {
      vector<vector<CCTK_REAL> > data; // [time level][var index]
    };
    // TODO: right now this is sized for the total number of groups
    vector<unique_ptr<ScalarGroupData> > scalargroupdata; // [group index]
  };
  GlobalData globaldata;

  struct LevelData {
    int level;
    // Empty MultiFab holding a cell-centred BoxArray for iterating
    // over grid functions.
    // TODO: Can we store the BoxArray directly?
    // The data in mfab0 is not valid. It's a dummy
    // variable to get the distribution mapping, the
    // grid size, ghost zones, etc. Only needed for
    // cctkGH.
    unique_ptr<MultiFab> mfab0;

    struct GroupData : public CommonGroupData {
      array<int, dim> indextype;
      array<int, dim> nghostzones;
      vector<array<int, dim> > parities;

      // each MultiFab has numvars components
      vector<unique_ptr<MultiFab> > mfab; // [time level]

      // flux register between this and the next coarser level
      unique_ptr<FluxRegister> freg;
      // associated flux group indices
      array<int, dim> fluxes; // [dir]
    };
    // TODO: right now this is sized for the total number of groups
    vector<unique_ptr<GroupData> > groupdata; // [group index]
  };
  vector<LevelData> leveldata; // [reflevel]
};

extern unique_ptr<GHExt> ghext;

Interpolater *get_interpolator(const array<int, dim> indextype);

typedef void apply_physbcs_t(const Box &, const FArrayBox &, int, int,
                             const Geometry &, CCTK_REAL, const Vector<BCRec> &,
                             int, int);
typedef PhysBCFunct<apply_physbcs_t *> CarpetXPhysBCFunct;
tuple<CarpetXPhysBCFunct, Vector<BCRec> >
get_boundaries(const GHExt::LevelData &leveldata,
               const GHExt::LevelData::GroupData &groupdata);

} // namespace CarpetX

#endif // #ifndef DRIVER_HXX
