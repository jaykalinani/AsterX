#include "driver.hxx"
#include "io.hxx"
#include "logo.hxx"
#include "prolongate_3d_rf2.hxx"
#include "schedule.hxx"
#include "timer.hxx"
#include "interp.hxx"

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <util_Table.h>

#include <AMReX.H>
#include <AMReX_BCRec.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PhysBCFunct.H>

#include <omp.h>
#include <mpi.h>

#include <array>
#include <cassert>
#include <cmath>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

namespace CarpetX {
using namespace std;

template <typename T> constexpr T ipow(T x, int n) {
  assert(n >= 0);
  T r = 1;
  while (n > 0) {
    if (n % 2)
      r *= x;
    n /= 2;
    x *= x;
  }
  return r;
}

// Global variables

int ghext_handle = -1;

amrex::AMReX *restrict pamrex = nullptr;
unique_ptr<GHExt> ghext;

// Registered functions

void *SetupGH(tFleshConfig *fc, int convLevel, cGH *cctkGH);
int InitGH(cGH *cctkGH);
int ScheduleTraverseGH(cGH *cctkGH, const char *where);

int MyProc(const cGH *cctkGH);
int nProcs(const cGH *cctkGH);
int Exit(cGH *cctkGH, int retval);
int Abort(cGH *cctkGH, int retval);
int Barrier(const cGH *cctkGHa);
int GroupDynamicData(const cGH *cctkGH, int gi, cGroupDynamicData *data);

// Aliased functions

extern "C" void CarpetX_CallScheduleGroup(void *cctkGH, const char *groupname);
extern "C" CCTK_INT CarpetX_GetDomainSpecification(
    CCTK_INT size, CCTK_REAL *restrict const physical_min,
    CCTK_REAL *restrict const physical_max,
    CCTK_REAL *restrict const interior_min,
    CCTK_REAL *restrict const interior_max,
    CCTK_REAL *restrict const exterior_min,
    CCTK_REAL *restrict const exterior_max, CCTK_REAL *restrict const spacing);
extern "C" CCTK_INT CarpetX_GetBoundarySizesAndTypes(
    const void *cctkGH_, CCTK_INT size, CCTK_INT *restrict const bndsize,
    CCTK_INT *restrict const is_ghostbnd, CCTK_INT *restrict const is_symbnd,
    CCTK_INT *restrict const is_physbnd);

// Local functions
void SetupGlobals();

////////////////////////////////////////////////////////////////////////////////

bool get_group_checkpoint_flag(const int gi) {
  int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  char buf[100];
  int iret = Util_TableGetString(tags, sizeof buf, buf, "checkpoint");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    return true;
  } else if (iret >= 0) {
    string str(buf);
    for (auto &c : str)
      c = tolower(c);
    if (str == "yes")
      return true;
    if (str == "no")
      return false;
    assert(0);
  } else {
    assert(0);
  }
}

bool get_group_restrict_flag(const int gi) {
  int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  char buf[100];
  int iret = Util_TableGetString(tags, sizeof buf, buf, "restrict");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    return true;
  } else if (iret >= 0) {
    string str(buf);
    for (auto &c : str)
      c = tolower(c);
    if (str == "yes")
      return true;
    if (str == "no")
      return false;
    assert(0);
  } else {
    assert(0);
  }
}

array<int, dim> get_group_indextype(const int gi) {
  assert(gi >= 0);
  int tags = CCTK_GroupCenteringTableI(gi);
  const char *tag_name = "centering";
  if(tags < 0) {
      tags = CCTK_GroupTagsTableI(gi);
      tag_name = "index";
  }
  assert(tags >= 0);
  array<CCTK_INT, dim> index;
  int iret = Util_TableGetIntArray(tags, dim, index.data(), tag_name);
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    index = {0, 0, 0}; // default: vertex-centred
  } else if (iret >= 0) {
    assert(iret == dim);
  } else {
    assert(0);
  }
  array<int, dim> indextype;
  for (int d = 0; d < dim; ++d)
    indextype[d] = index[d];
  return indextype;
}

array<int, dim> get_group_fluxes(const int gi) {
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  vector<char> fluxes_buf(1000);
  const int iret =
      Util_TableGetString(tags, fluxes_buf.size(), fluxes_buf.data(), "fluxes");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    fluxes_buf[0] = '\0'; // default: empty (no fluxes)
  } else if (iret >= 0) {
    // do nothing
  } else {
    assert(0);
  }

  const string str(fluxes_buf.data());
  vector<string> strs;
  size_t end = 0;
  while (end < str.size()) {
    const size_t begin = str.find_first_not_of(' ', end);
    if (begin == string::npos)
      break;
    end = str.find(' ', begin);
    strs.push_back(str.substr(begin, end - begin));
  }

  array<int, dim> fluxes;
  fluxes.fill(-1);
  if (strs.empty())
    return fluxes; // No fluxes specified

  assert(strs.size() == dim); // Check number of fluxes
  for (int d = 0; d < dim; ++d) {
    auto str1 = strs[d];
    if (str1.find(':') == string::npos) {
      const char *const impl = CCTK_GroupImplementationI(gi);
      str1 = string(impl) + "::" + str1;
    }
    int gi1 = CCTK_GroupIndex(str1.c_str());
    assert(gi1 >= 0); // Check fluxes are valid groups
    fluxes[d] = gi1;
  }

  for (int d = 0; d < dim; ++d)
    for (int d1 = d + 1; d1 < dim; ++d1)
      assert(fluxes[d] != fluxes[d1]); // Check groups are all different
  return fluxes;
}

array<int, dim> get_group_nghostzones(const int gi) {
  DECLARE_CCTK_PARAMETERS;
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  array<CCTK_INT, dim> nghostzones;
  int iret =
      Util_TableGetIntArray(tags, dim, nghostzones.data(), "nghostzones");
  if (iret == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    // default: use driver parameter
    nghostzones = {ghost_size >= 0 ? ghost_size : ghost_size_x,
                   ghost_size >= 0 ? ghost_size : ghost_size_y,
                   ghost_size >= 0 ? ghost_size : ghost_size_z};
  } else if (iret >= 0) {
    assert(iret == dim);
  } else {
    assert(0);
  }
  return nghostzones;
}

vector<array<int, dim> > get_group_parities(const int gi) {
  DECLARE_CCTK_PARAMETERS;
  assert(gi >= 0);
  const int tags = CCTK_GroupTagsTableI(gi);
  assert(tags >= 0);
  const int nelems = Util_TableGetIntArray(tags, 0, nullptr, "parities");
  if (nelems == UTIL_ERROR_TABLE_NO_SUCH_KEY) {
    // unset (will use index type)
    return {};
  } else if (nelems >= 0) {
    // do nothing
  } else {
    assert(0);
  }
  vector<CCTK_INT> parities1(nelems);
  const int iret =
      Util_TableGetIntArray(tags, nelems, parities1.data(), "parities");
  assert(iret == nelems);
  assert(nelems % dim == 0);
  vector<array<int, dim> > parities(nelems / dim);
  for (size_t n = 0; n < parities.size(); ++n)
    for (int d = 0; d < dim; ++d)
      parities.at(n)[d] = parities1.at(dim * n + d);
  return parities;
}

amrex::Interpolater *get_interpolator(const array<int, dim> indextype) {
  DECLARE_CCTK_PARAMETERS;

  enum class interp_t { unset, interpolate, conservative, ddf };
  static interp_t interp = [&]() {
    if (CCTK_EQUALS(prolongation_type, "interpolate")) {
      return interp_t::interpolate;
    } else if (CCTK_EQUALS(prolongation_type, "conservative")) {
      return interp_t::conservative;
    } else if (CCTK_EQUALS(prolongation_type, "ddf")) {
      return interp_t::ddf;
    } else {
      assert(0);
    }
  }();

  switch (interp) {
  case interp_t::interpolate:

    switch ((indextype[0] << 2) | (indextype[1] << 1) | (indextype[2] << 0)) {

    case 0b000:
      switch (prolongation_order) {
      case 1:
        return &prolongate_3d_rf2_c000_o1;
      case 3:
        return &prolongate_3d_rf2_c000_o3;
      }
      break;

    case 0b001:
      switch (prolongation_order) {
      case 1:
        return &prolongate_3d_rf2_c001_o1;
      case 3:
        return &prolongate_3d_rf2_c001_o3;
      }
      break;

    case 0b010:
      switch (prolongation_order) {
      case 1:
        return &prolongate_3d_rf2_c010_o1;
      case 3:
        return &prolongate_3d_rf2_c010_o3;
      }
      break;

    case 0b011:
      switch (prolongation_order) {
      case 1:
        return &prolongate_3d_rf2_c011_o1;
      case 3:
        return &prolongate_3d_rf2_c011_o3;
      }
      break;

    case 0b100:
      switch (prolongation_order) {
      case 1:
        return &prolongate_3d_rf2_c100_o1;
      case 3:
        return &prolongate_3d_rf2_c100_o3;
      }
      break;

    case 0b101:
      switch (prolongation_order) {
      case 1:
        return &prolongate_3d_rf2_c101_o1;
      case 3:
        return &prolongate_3d_rf2_c101_o3;
      }
      break;

    case 0b110:
      switch (prolongation_order) {
      case 1:
        return &prolongate_3d_rf2_c110_o1;
      case 3:
        return &prolongate_3d_rf2_c110_o3;
      }
      break;

    case 0b111:
      switch (prolongation_order) {
      case 1:
        return &prolongate_3d_rf2_c111_o1;
      case 3:
        return &prolongate_3d_rf2_c111_o3;
      }
      break;
    }
    break;

  case interp_t::conservative:

    switch ((indextype[0] << 2) | (indextype[1] << 1) | (indextype[2] << 0)) {

    case 0b000:
      switch (prolongation_order) {
      case 1:
        return &prolongate_cons_3d_rf2_c000_o0;
      }
      break;

    case 0b001:
      switch (prolongation_order) {
      case 1:
        return &prolongate_cons_3d_rf2_c001_o0;
      }
      break;

    case 0b010:
      switch (prolongation_order) {
      case 1:
        return &prolongate_cons_3d_rf2_c010_o0;
      }
      break;

    case 0b011:
      switch (prolongation_order) {
      case 1:
        return &prolongate_cons_3d_rf2_c011_o0;
      }
      break;

    case 0b100:
      switch (prolongation_order) {
      case 1:
        return &prolongate_cons_3d_rf2_c100_o0;
      }
      break;

    case 0b101:
      switch (prolongation_order) {
      case 1:
        return &prolongate_cons_3d_rf2_c101_o0;
      }
      break;

    case 0b110:
      switch (prolongation_order) {
      case 1:
        return &prolongate_cons_3d_rf2_c110_o0;
      }
      break;

    case 0b111:
      switch (prolongation_order) {
      case 1:
        return &prolongate_cons_3d_rf2_c111_o0;
      }
      break;
    }
    break;

  case interp_t::ddf:

    switch (prolongation_order) {

    case 1:
      switch ((indextype[0] << 2) | (indextype[1] << 1) | (indextype[2] << 0)) {
      case 0b000:
        return &prolongate_ddf_3d_rf2_c000_o1;
      case 0b001:
        return &prolongate_ddf_3d_rf2_c001_o1;
      case 0b010:
        return &prolongate_ddf_3d_rf2_c010_o1;
      case 0b011:
        return &prolongate_ddf_3d_rf2_c011_o1;
      case 0b100:
        return &prolongate_ddf_3d_rf2_c100_o1;
      case 0b101:
        return &prolongate_ddf_3d_rf2_c101_o1;
      case 0b110:
        return &prolongate_ddf_3d_rf2_c110_o1;
      case 0b111:
        return &prolongate_ddf_3d_rf2_c111_o1;
      }
      break;

    case 3:
      switch ((indextype[0] << 2) | (indextype[1] << 1) | (indextype[2] << 0)) {
      case 0b000:
        return &prolongate_ddf_3d_rf2_c000_o3;
      case 0b001:
        return &prolongate_ddf_3d_rf2_c001_o3;
      case 0b010:
        return &prolongate_ddf_3d_rf2_c010_o3;
      case 0b011:
        return &prolongate_ddf_3d_rf2_c011_o3;
      case 0b100:
        return &prolongate_ddf_3d_rf2_c100_o3;
      case 0b101:
        return &prolongate_ddf_3d_rf2_c101_o3;
      case 0b110:
        return &prolongate_ddf_3d_rf2_c110_o3;
      case 0b111:
        return &prolongate_ddf_3d_rf2_c111_o3;
      }
      break;

    case 5:
      switch ((indextype[0] << 2) | (indextype[1] << 1) | (indextype[2] << 0)) {
      case 0b000:
        return &prolongate_ddf_3d_rf2_c000_o5;
      case 0b001:
        return &prolongate_ddf_3d_rf2_c001_o5;
      case 0b010:
        return &prolongate_ddf_3d_rf2_c010_o5;
      case 0b011:
        return &prolongate_ddf_3d_rf2_c011_o5;
      case 0b100:
        return &prolongate_ddf_3d_rf2_c100_o5;
      case 0b101:
        return &prolongate_ddf_3d_rf2_c101_o5;
      case 0b110:
        return &prolongate_ddf_3d_rf2_c110_o5;
      case 0b111:
        return &prolongate_ddf_3d_rf2_c111_o5;
      }
      break;
    }
    break;

  case interp_t::unset:
    // do nothing; errors are handled below
    break;

  } // switch prolongation_type

  CCTK_VERROR("Unsupported combination of prolongation_type \"%s\", "
              "prolongation order %d, and index type [%d,%d,%d]",
              prolongation_type, prolongation_order, indextype[0], indextype[1],
              indextype[2]);
  assert(0);
}

void apply_physbcs(const amrex::Box &, const amrex::FArrayBox &, int, int,
                   const amrex::Geometry &, CCTK_REAL,
                   const amrex::Vector<amrex::BCRec> &, int,
                   int) { // do nothing
}

tuple<CarpetXPhysBCFunct, amrex::Vector<amrex::BCRec> >
get_boundaries(const GHExt::LevelData::GroupData &groupdata) {
  DECLARE_CCTK_PARAMETERS;

  // TODO: It seems that AMReX now also has `RB90`, `RB180`, and
  // `PolarB` boundary conditions. Make these available as well.
  const array<array<bool, 3>, 2> is_periodic{{
      {{
          periodic || periodic_x,
          periodic || periodic_y,
          periodic || periodic_z,
      }},
      {{
          periodic || periodic_x,
          periodic || periodic_y,
          periodic || periodic_z,
      }},
  }};
  const array<array<bool, 3>, 2> is_reflect{{
      {{
          !!reflection_x,
          !!reflection_y,
          !!reflection_z,
      }},
      {{
          !!reflection_upper_x,
          !!reflection_upper_y,
          !!reflection_upper_z,
      }},
  }};
  const auto makebc = [&](const int vi, const int dir, const int face) {
    assert(dir >= 0 && dir < dim);
    assert(face >= 0 && face < 2);
    if (is_periodic[face][dir])
      return amrex::BCType::int_dir;
    if (is_reflect[face][dir])
      return groupdata.parities.at(vi)[dir] > 0 ? amrex::BCType::reflect_even
                                                : amrex::BCType::reflect_odd;
    return amrex::BCType::ext_dir;
  };

  amrex::Vector<amrex::BCRec> bcs(groupdata.numvars);
  for (int vi = 0; vi < groupdata.numvars; ++vi)
    bcs.at(vi) =
        amrex::BCRec(makebc(vi, 0, 0), makebc(vi, 1, 0), makebc(vi, 2, 0),
                     makebc(vi, 0, 1), makebc(vi, 1, 1), makebc(vi, 2, 1));
  CarpetXPhysBCFunct physbc(ghext->amrcore->Geom(groupdata.level), bcs,
                            apply_physbcs);

  return {move(physbc), move(bcs)};
}

////////////////////////////////////////////////////////////////////////////////

// amrex::AmrCore functions

CactusAmrCore::CactusAmrCore() {}
CactusAmrCore::CactusAmrCore(const amrex::RealBox *rb, int max_level_in,
                             const amrex::Vector<int> &n_cell_in, int coord,
                             amrex::Vector<amrex::IntVect> ref_ratios,
                             const int *is_per)
    : amrex::AmrCore(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per) {
  SetupGlobals();
}

CactusAmrCore::CactusAmrCore(const amrex::RealBox &rb, int max_level_in,
                             const amrex::Vector<int> &n_cell_in, int coord,
                             amrex::Vector<amrex::IntVect> const &ref_ratios,
                             amrex::Array<int, AMREX_SPACEDIM> const &is_per)
    : amrex::AmrCore(rb, max_level_in, n_cell_in, coord, ref_ratios, is_per) {
  SetupGlobals();
}

CactusAmrCore::~CactusAmrCore() {}

void CactusAmrCore::ErrorEst(const int level, amrex::TagBoxArray &tags,
                             const amrex::Real time, const int ngrow) {
  DECLARE_CCTK_PARAMETERS;

  // Don't regrid before Cactus is ready to
  if (level >= int(ghext->leveldata.size()))
    return;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("ErrorEst level %d", level);

  const int gi = CCTK_GroupIndex("CarpetX::regrid_error");
  assert(gi >= 0);
  const int vi = 0;
  const int tl = 0;

  auto &restrict leveldata = ghext->leveldata.at(level);
  auto &restrict groupdata = *leveldata.groupdata.at(gi);
  // Ensure the error estimate has been set
  error_if_invalid(groupdata, vi, tl, make_valid_int(),
                   []() { return "ErrorEst"; });
  auto mfitinfo = amrex::MFItInfo().SetDynamic(true).EnableTiling();
#pragma omp parallel
  for (amrex::MFIter mfi(*leveldata.fab, mfitinfo); mfi.isValid(); ++mfi) {
    GridPtrDesc1 grid(groupdata, mfi);

    const amrex::Array4<const CCTK_REAL> &err_array4 =
        groupdata.mfab.at(tl)->array(mfi);
    const GF3D1<const CCTK_REAL> &err_ = grid.gf3d(err_array4, vi);
    const amrex::Array4<char> &tags_array4 = tags.array(mfi);

    grid.loop_idx(
        where_t::interior, groupdata.indextype, [&](const Loop::PointDesc &p) {
          tags_array4(grid.cactus_offset.x + p.i, grid.cactus_offset.y + p.j,
                      grid.cactus_offset.z + p.k) =
              err_(p.I) >= regrid_error_threshold ? amrex::TagBox::SET
                                                  : amrex::TagBox::CLEAR;
        });
    // Do not set the boundary; AMReX's error grid function might have
    // a different number of ghost zones, and these ghost zones are
    // probably unused anyway.
    if (false) {
      grid.loop_idx(
          where_t::boundary, groupdata.indextype,
          [&](const Loop::PointDesc &p) {
            tags_array4(grid.cactus_offset.x + p.i, grid.cactus_offset.y + p.j,
                        grid.cactus_offset.z + p.k) = amrex::TagBox::CLEAR;
          });
    }
  }
}

void SetupGlobals() {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("SetupGlobals");

  GHExt::GlobalData &globaldata = ghext->globaldata;

  const int numgroups = CCTK_NumGroups();
  globaldata.arraygroupdata.resize(numgroups);
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    // grid functions live on levels
    if (group.grouptype == CCTK_GF)
      continue;
    assert(group.grouptype == CCTK_ARRAY || group.grouptype == CCTK_SCALAR);
    assert(group.vartype == CCTK_VARIABLE_REAL);
    assert(group.disttype == CCTK_DISTRIB_CONSTANT);
    assert(group.dim >= 0);

    globaldata.arraygroupdata.at(gi) =
        make_unique<GHExt::GlobalData::ArrayGroupData>();
    GHExt::GlobalData::ArrayGroupData &arraygroupdata =
        *globaldata.arraygroupdata.at(gi);
    arraygroupdata.groupindex = gi;
    arraygroupdata.firstvarindex = CCTK_FirstVarIndexI(gi);
    arraygroupdata.numvars = group.numvars;
    arraygroupdata.do_checkpoint = get_group_checkpoint_flag(gi);
    arraygroupdata.do_restrict = get_group_restrict_flag(gi);

    CCTK_INT const *const *const sz = CCTK_GroupSizesI(gi);
    arraygroupdata.array_size = 1;
    for (int d = 0; d < group.dim; ++d) {
      arraygroupdata.array_size = arraygroupdata.array_size * *sz[d];
    }

    // Set up dynamic data
    arraygroupdata.dimension = group.dim;
    arraygroupdata.activetimelevels = 1;
    for (int d = 0; d < group.dim; ++d) {
      arraygroupdata.lsh[d] = *sz[d];
      arraygroupdata.ash[d] = *sz[d];
      arraygroupdata.gsh[d] = *sz[d];
      arraygroupdata.nghostzones[d] = 0;
      arraygroupdata.lbnd[d] = 0;
      arraygroupdata.ubnd[d] = *sz[d] - 1;
      arraygroupdata.bbox[2 * d] = arraygroupdata.bbox[2 * d + 1] = 1;
    }

    // Allocate data
    arraygroupdata.data.resize(group.numtimelevels);
    arraygroupdata.valid.resize(group.numtimelevels);
    for (int tl = 0; tl < int(arraygroupdata.data.size()); ++tl) {
      arraygroupdata.data.at(tl).resize(arraygroupdata.numvars *
                                        arraygroupdata.array_size);
      why_valid_t why([]() { return "SetupGlobals"; });
      arraygroupdata.valid.at(tl).resize(arraygroupdata.numvars, why);
      for (int vi = 0; vi < arraygroupdata.numvars; ++vi) {
        // TODO: decide that valid_bnd == false always and rely on
        // initialization magic?
        valid_t valid;
        valid.valid_int = false;
        valid.valid_outer = true;
        valid.valid_ghosts = true;
        arraygroupdata.valid.at(tl).at(vi).set(valid,
                                               []() { return "SetupGlobals"; });

        // TODO: make poison_invalid and check_invalid virtual members of
        // CommonGroupData
        poison_invalid(arraygroupdata, vi, tl);
        check_valid(arraygroupdata, vi, tl, []() { return "SetupGlobals"; });
      }
    }
  }
}

void SetupLevel(const int level, const amrex::BoxArray &ba,
                const amrex::DistributionMapping &dm,
                const function<string()> &why) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("SetupLevel level %d", level);

  assert(level == int(ghext->leveldata.size()));

  ghext->leveldata.resize(level + 1);
  GHExt::LevelData &leveldata = ghext->leveldata.at(level);
  leveldata.level = level;

  const int timereffact = use_subcycling_wip ? 2 : 1;
  if (level == 0) {
    // We are creating the coarsest level
    leveldata.is_subcycling_level = false; // unused
    leveldata.iteration = 0;
    leveldata.delta_iteration = 1;
  } else {
    // We are creating a new refined level
    auto &coarseleveldata = ghext->leveldata.at(level - 1);
    leveldata.is_subcycling_level = use_subcycling_wip;
    leveldata.iteration = coarseleveldata.iteration;
    leveldata.delta_iteration = coarseleveldata.delta_iteration / timereffact;
  }

  const amrex::IntVect nghostzones = {
      ghost_size >= 0 ? ghost_size : ghost_size_x,
      ghost_size >= 0 ? ghost_size : ghost_size_y,
      ghost_size >= 0 ? ghost_size : ghost_size_z};
  leveldata.fab = make_unique<amrex::FabArrayBase>(ba, dm, 1, nghostzones);
  assert(ba.ixType() == amrex::IndexType(amrex::IndexType::CELL,
                                         amrex::IndexType::CELL,
                                         amrex::IndexType::CELL));

  const auto &geom = ghext->amrcore->Geom(leveldata.level);

  const int numgroups = CCTK_NumGroups();
  leveldata.groupdata.resize(numgroups);
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    // only grid functions live on levels (and the grid)
    if (group.grouptype != CCTK_GF)
      continue;

    assert(group.grouptype == CCTK_GF);
    assert(group.vartype == CCTK_VARIABLE_REAL);
    assert(group.disttype == CCTK_DISTRIB_DEFAULT);
    assert(group.dim == dim);

    leveldata.groupdata.at(gi) = make_unique<GHExt::LevelData::GroupData>();
    GHExt::LevelData::GroupData &groupdata = *leveldata.groupdata.at(gi);
    groupdata.level = leveldata.level;
    groupdata.groupindex = gi;
    groupdata.firstvarindex = CCTK_FirstVarIndexI(gi);
    groupdata.numvars = group.numvars;
    groupdata.do_checkpoint = get_group_checkpoint_flag(gi);
    groupdata.do_restrict = get_group_restrict_flag(gi);
    groupdata.indextype = get_group_indextype(gi);
    groupdata.nghostzones = get_group_nghostzones(gi);
    groupdata.parities = get_group_parities(gi);
    if (groupdata.parities.empty()) {
      array<int, dim> parity;
      for (int d = 0; d < dim; ++d)
        parity[d] = groupdata.indextype[d] == 0 ? +1 : -1;
      groupdata.parities.resize(groupdata.numvars, parity);
    }
    assert(int(groupdata.parities.size()) == groupdata.numvars);
    for (int vi = 0; vi < groupdata.numvars; ++vi)
      for (int d = 0; d < dim; ++d)
        assert(abs(groupdata.parities.at(vi)[d]) == 1);

    // Periodic boundaries require (num interior points) >= (num ghost points)
    for (int d = 0; d < dim; ++d)
      if (geom.isPeriodic(d))
        assert(geom.Domain().length(d) >= groupdata.nghostzones[d]);

    // Allocate grid hierarchies
    // Note: CarpetX and AMReX use opposite constants for cell/vertex
    // data. AMReX uses (AMReX_IndexType.H): enum CellIndex { CELL =
    // 0, NODE = 1 }, while CarpetX uses indextype == 0 for vertex
    // (node) and indextype == 1 for cell data.
    const amrex::BoxArray &gba = convert(
        ba, amrex::IndexType(groupdata.indextype[0] ? amrex::IndexType::CELL
                                                    : amrex::IndexType::NODE,
                             groupdata.indextype[1] ? amrex::IndexType::CELL
                                                    : amrex::IndexType::NODE,
                             groupdata.indextype[2] ? amrex::IndexType::CELL
                                                    : amrex::IndexType::NODE));
    groupdata.mfab.resize(group.numtimelevels);
    groupdata.valid.resize(group.numtimelevels);
    for (int tl = 0; tl < int(groupdata.mfab.size()); ++tl) {
      groupdata.mfab.at(tl) = make_unique<amrex::MultiFab>(
          gba, dm, groupdata.numvars, amrex::IntVect(groupdata.nghostzones));
      groupdata.valid.at(tl).resize(groupdata.numvars, why_valid_t(why));
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        poison_invalid(groupdata, vi, tl);
    }

    if (level > 0) {
      array<int, dim> fluxes = get_group_fluxes(groupdata.groupindex);
      const bool have_fluxes = fluxes[0] >= 0;
      if (have_fluxes) {
        assert((groupdata.indextype == array<int, dim>{1, 1, 1}));
        groupdata.freg = make_unique<amrex::FluxRegister>(
            gba, dm, ghext->amrcore->refRatio(level - 1), level,
            groupdata.numvars);
        groupdata.fluxes = fluxes;
      } else {
        groupdata.fluxes.fill(-1);
      }
    }
  }

  // Check flux register consistency
  for (int gi = 0; gi < numgroups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    /* only grid functions live on levels (and the grid) */
    if (group.grouptype != CCTK_GF)
      continue;
    const auto &groupdata = *leveldata.groupdata.at(gi);
    if (groupdata.freg) {
      for (int d = 0; d < dim; ++d) {
        assert(groupdata.fluxes[d] != groupdata.groupindex);
        const auto &flux_groupdata =
            *leveldata.groupdata.at(groupdata.fluxes[d]);
        array<int, dim> flux_indextype{1, 1, 1};
        flux_indextype[d] = 0;
        assert(flux_groupdata.indextype == flux_indextype);
        assert(flux_groupdata.numvars == groupdata.numvars);
      }
    }
  }
}

void CactusAmrCore::MakeNewLevelFromScratch(
    int level, amrex::Real time, const amrex::BoxArray &ba,
    const amrex::DistributionMapping &dm) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("MakeNewLevelFromScratch level %d", level);

  SetupLevel(level, ba, dm, []() { return "MakeNewLevelFromScratch"; });

  cGH *const cctkGH = get_global_cctkGH();
  if (cctkGH) {
    setup_cctkGHs(cctkGH);
    assert(!active_levels);
    active_levels = make_optional<active_levels_t>(level, level + 1);
    CCTK_Traverse(cctkGH, "CCTK_BASEGRID");
    // CCTK_Traverse(cctkGH, "CCTK_POSTREGRID");
    active_levels = optional<active_levels_t>();
  }
}

void CactusAmrCore::MakeNewLevelFromCoarse(
    const int level, const amrex::Real time, const amrex::BoxArray &ba,
    const amrex::DistributionMapping &dm) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("MakeNewLevelFromCoarse level %d", level);

  assert(!use_subcycling_wip);
  assert(level > 0);

  SetupLevel(level, ba, dm, []() { return "MakeNewLevelFromCoarse"; });

  // Prolongate
  assert(!use_subcycling_wip);
  auto &leveldata = ghext->leveldata.at(level);
  auto &coarseleveldata = ghext->leveldata.at(level - 1);
  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata = *leveldata.groupdata.at(gi);
    auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
    assert(coarsegroupdata.numvars == groupdata.numvars);
    amrex::Interpolater *const interpolator =
        get_interpolator(groupdata.indextype);

    const amrex::IntVect reffact{2, 2, 2};

    auto physbc_bcs = get_boundaries(groupdata);
    CarpetXPhysBCFunct &physbc = get<0>(physbc_bcs);
    const amrex::Vector<amrex::BCRec> &bcs = get<1>(physbc_bcs);

    const int ntls = groupdata.mfab.size();
    // We only prolongate the state vector. And if there is more than
    // one time level, then we don't prolongate the oldest.
    const int prolongate_tl =
        groupdata.do_checkpoint ? (ntls > 1 ? ntls - 1 : ntls) : 0;

    groupdata.valid.resize(ntls);
    for (int tl = 0; tl < ntls; ++tl) {
      why_valid_t why([]() { return "MakeNewLevelFromCoarse"; });
      groupdata.valid.at(tl).resize(groupdata.numvars, why);
      for (int vi = 0; vi < groupdata.numvars; ++vi)
        groupdata.valid.at(tl).at(vi).set(valid_t(false), []() {
          return "MakeNewLevelFromCoarse: not prolongated because variable is "
                 "not evolved";
        });

      if (tl < prolongate_tl) {
        // Expect coarse grid data to be valid
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          error_if_invalid(coarsegroupdata, vi, tl, make_valid_all(), []() {
            return "MakeNewLevelFromCoarse before prolongation";
          });
          check_valid(coarsegroupdata, vi, tl, []() {
            return "MakeNewLevelFromCoarse before prolongation";
          });
        }
        InterpFromCoarseLevel(
            *groupdata.mfab.at(tl), 0.0, *coarsegroupdata.mfab.at(tl), 0, 0,
            groupdata.numvars, ghext->amrcore->Geom(level - 1),
            ghext->amrcore->Geom(level), physbc, 0, physbc, 0, reffact,
            interpolator, bcs, 0);
        for (int vi = 0; vi < groupdata.numvars; ++vi)
          groupdata.valid.at(tl).at(vi).set(make_valid_int(), []() {
            return "MakeNewLevelFromCoarse after prolongation";
          });
      }

      for (int vi = 0; vi < groupdata.numvars; ++vi) {
        // Already poisoned by SetupLevel
        check_valid(groupdata, vi, tl, []() {
          return "MakeNewLevelFromCoarse after prolongation";
        });
      }
    } // for tl

  } // for gi

  cGH *const cctkGH = get_global_cctkGH();
  if (cctkGH) {
    setup_cctkGHs(cctkGH);
    assert(!active_levels);
    active_levels = make_optional<active_levels_t>(level, level + 1);
    CCTK_Traverse(cctkGH, "CCTK_BASEGRID");
    CCTK_Traverse(cctkGH, "CCTK_POSTREGRID");
    active_levels = optional<active_levels_t>();
  }
}

void CactusAmrCore::RemakeLevel(const int level, const amrex::Real time,
                                const amrex::BoxArray &ba,
                                const amrex::DistributionMapping &dm) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("RemakeLevel level %d", level);

  // Copy or prolongate
  auto &leveldata = ghext->leveldata.at(level);
  assert(!use_subcycling_wip);
  assert(leveldata.level > 0);
  auto &coarseleveldata = ghext->leveldata.at(level - 1);

  const amrex::IntVect nghostzones = {
      ghost_size >= 0 ? ghost_size : ghost_size_x,
      ghost_size >= 0 ? ghost_size : ghost_size_y,
      ghost_size >= 0 ? ghost_size : ghost_size_z};
  leveldata.fab = make_unique<amrex::FabArrayBase>(ba, dm, 1, nghostzones);
  assert(ba.ixType() == amrex::IndexType(amrex::IndexType::CELL,
                                         amrex::IndexType::CELL,
                                         amrex::IndexType::CELL));

  // We assume that this level is at the same time as the next coarser level
  assert(leveldata.iteration == coarseleveldata.iteration);

  const int num_groups = CCTK_NumGroups();
  for (int gi = 0; gi < num_groups; ++gi) {
    cGroup group;
    int ierr = CCTK_GroupData(gi, &group);
    assert(!ierr);

    if (group.grouptype != CCTK_GF)
      continue;

    auto &restrict groupdata = *leveldata.groupdata.at(gi);
    auto &restrict coarsegroupdata = *coarseleveldata.groupdata.at(gi);
    assert(coarsegroupdata.numvars == groupdata.numvars);
    amrex::Interpolater *const interpolator =
        get_interpolator(groupdata.indextype);

    const amrex::IntVect reffact{2, 2, 2};

    auto physbc_bcs = get_boundaries(groupdata);
    CarpetXPhysBCFunct &physbc = get<0>(physbc_bcs);
    const amrex::Vector<amrex::BCRec> &bcs = get<1>(physbc_bcs);

    const amrex::BoxArray &gba = convert(
        ba, amrex::IndexType(groupdata.indextype[0] ? amrex::IndexType::CELL
                                                    : amrex::IndexType::NODE,
                             groupdata.indextype[1] ? amrex::IndexType::CELL
                                                    : amrex::IndexType::NODE,
                             groupdata.indextype[2] ? amrex::IndexType::CELL
                                                    : amrex::IndexType::NODE));

    const int ntls = groupdata.mfab.size();
    // We only prolongate the state vector. And if there is more than
    // one time level, then we don't prolongate the oldest.
    const int prolongate_tl =
        groupdata.do_checkpoint ? (ntls > 1 ? ntls - 1 : ntls) : 0;

    for (int tl = 0; tl < ntls; ++tl) {
      auto mfab = make_unique<amrex::MultiFab>(
          gba, dm, groupdata.numvars, amrex::IntVect(groupdata.nghostzones));
      vector<why_valid_t> valid(groupdata.numvars, why_valid_t([]() {
                                  return "RemakeLevel: not prolongated/copied "
                                         "because variable is not evolved";
                                }));
      // We cannot poison the grid variables because the grid
      // structure in `leveldata.fab` already changed, but the new
      // variables in `mfab` have not yet been assigned to the driver
      // data structure.
      // for (int vi = 0; vi < groupdata.numvars; ++vi)
      //   poison_invalid(groupdata, vi, tl);

      if (tl < prolongate_tl) {
        for (int vi = 0; vi < groupdata.numvars; ++vi) {
          error_if_invalid(coarsegroupdata, vi, tl, make_valid_all(),
                           []() { return "RemakeLevel before prolongation"; });
          error_if_invalid(groupdata, vi, tl, make_valid_all(),
                           []() { return "RemakeLevel before prolongation"; });
          check_valid(coarsegroupdata, vi, tl,
                      []() { return "RemakeLevel before prolongation"; });
          // We cannot call this function since it would try to
          // traverse the old grid function with the new grid
          // structure.
          // TODO: Use explicit loop instead (see above)
          // check_valid(leveldata, groupdata, vi, tl);
        }

        // Copy from same level and/or prolongate from next coarser level
        FillPatchTwoLevels(*mfab, 0.0, {&*coarsegroupdata.mfab.at(tl)}, {0.0},
                           {&*groupdata.mfab.at(tl)}, {0.0}, 0, 0,
                           groupdata.numvars, ghext->amrcore->Geom(level - 1),
                           ghext->amrcore->Geom(level), physbc, 0, physbc, 0,
                           reffact, interpolator, bcs, 0);

        for (int vi = 0; vi < groupdata.numvars; ++vi)
          valid.at(vi) = why_valid_t(make_valid_int(), []() {
            return "RemakeLevel after prolongation";
          });
      }

      groupdata.mfab.at(tl) = move(mfab);
      groupdata.valid.at(tl) = move(valid);
      if (groupdata.freg)
        groupdata.freg = make_unique<amrex::FluxRegister>(
            gba, dm, ghext->amrcore->refRatio(level - 1), level,
            groupdata.numvars);

      for (int vi = 0; vi < groupdata.numvars; ++vi) {
        // Already poisoned above
        check_valid(groupdata, vi, tl,
                    []() { return "RemakeLevel after prolongation"; });
      }
    } // for tl

  } // for gi

  cGH *const cctkGH = get_global_cctkGH();
  if (cctkGH) {
    setup_cctkGHs(cctkGH);
    assert(!active_levels);
    active_levels = make_optional<active_levels_t>(level, level + 1);
    CCTK_Traverse(cctkGH, "CCTK_BASEGRID");
    CCTK_Traverse(cctkGH, "CCTK_POSTREGRID");
    active_levels = optional<active_levels_t>();
  }
}

void CactusAmrCore::ClearLevel(const int level) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("ClearLevel level %d", level);

  // note: several levels can be removed at once
  ghext->leveldata.resize(level);
}

////////////////////////////////////////////////////////////////////////////////

template <typename T, size_t N> inline vector<T> seq(const array<T, N> &v) {
  vector<T> r;
  for (const auto &x : v)
    r.push_back(x);
  return r;
}

template <typename T, size_t N>
inline vector<vector<T> > seqs(const vector<array<T, N> > &v) {
  vector<vector<T> > r;
  for (const auto &x : v)
    r.push_back(seq(x));
  return r;
}

} // namespace CarpetX

namespace std {
template <typename T, size_t N>
YAML::Emitter &operator<<(YAML::Emitter &yaml, const array<T, N> &arr) {
  yaml << YAML::Flow << YAML::BeginSeq;
  for (const auto &elt : arr)
    yaml << elt;
  yaml << YAML::EndSeq;
  return yaml;
}
} // namespace std

namespace amrex {

YAML::Emitter &operator<<(YAML::Emitter &yaml, const amrex::IntVect &iv) {
  yaml << YAML::Flow << YAML::BeginSeq;
  for (int d = 0; d < AMREX_SPACEDIM; ++d)
    yaml << iv[d];
  yaml << YAML::EndSeq;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const amrex::Box &box) {
  yaml << YAML::LocalTag("box-1.0.0");
  yaml << YAML::Flow << YAML::BeginMap;
  yaml << YAML::Key << "small" << YAML::Value << box.smallEnd();
  yaml << YAML::Key << "big" << YAML::Value << box.bigEnd();
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const amrex::RealBox &rbox) {
  yaml << YAML::LocalTag("realbox-1.0.0");
  yaml << YAML::Flow << YAML::BeginMap;
  yaml << YAML::Key << "xlo" << YAML::Value
       << std::vector<amrex::Real>(rbox.lo(), rbox.lo() + AMREX_SPACEDIM);
  yaml << YAML::Key << "xhi" << YAML::Value
       << std::vector<amrex::Real>(rbox.hi(), rbox.hi() + AMREX_SPACEDIM);
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const amrex::BoxArray &ba) {
  yaml << YAML::LocalTag("boxarray-1.0.0");
  yaml << YAML::BeginSeq;
  for (int n = 0; n < ba.size(); ++n)
    yaml << ba[n];
  yaml << YAML::EndSeq;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const amrex::Geometry &geom) {
  yaml << YAML::LocalTag("geometry-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "prob_domain" << YAML::Value << geom.ProbDomain();
  yaml << YAML::Key << "domain" << YAML::Value << geom.Domain();
  yaml << YAML::Key << "is_periodic" << YAML::Value << geom.isPeriodic();
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml,
                          const amrex::DistributionMapping &dm) {
  yaml << YAML::LocalTag("distributionmapping-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "processorMap" << YAML::Value << YAML::Flow
       << dm.ProcessorMap();
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const amrex::FabArrayBase &fab) {
  yaml << YAML::LocalTag("fabarraybase-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "ixType" << YAML::Value << fab.ixType().toIntVect();
  yaml << YAML::Key << "nGrowVect" << YAML::Value << fab.nGrowVect();
  yaml << YAML::Key << "boxArray" << YAML::Value << fab.boxArray();
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const amrex::AmrCore &amrcore) {
  yaml << YAML::LocalTag("amrcore-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "maxLevel" << YAML::Value << amrcore.maxLevel();
  yaml << YAML::Key << "finestLevel" << YAML::Value << amrcore.finestLevel();
  yaml << YAML::Key << "geometry" << YAML::Value << amrcore.Geom();
  yaml << YAML::Key << "distributionMapping" << YAML::Value
       << amrcore.DistributionMap();
  yaml << YAML::Key << "boxArray" << YAML::Value << amrcore.boxArray();
  yaml << YAML::EndMap;
  return yaml;
}

} // namespace amrex
namespace CarpetX {

YAML::Emitter &operator<<(YAML::Emitter &yaml,
                          const GHExt::CommonGroupData &commongroupdata) {
  yaml << YAML::LocalTag("commongroupdata-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "groupname" << YAML::Value
       << CCTK_FullGroupName(commongroupdata.groupindex);
  yaml << YAML::Key << "numvars" << YAML::Value << commongroupdata.numvars;
  yaml << YAML::Key << "varnames" << YAML::Value << YAML::Flow
       << YAML::BeginSeq;
  for (int vi = 0; vi < commongroupdata.numvars; ++vi)
    yaml << CCTK_VarName(commongroupdata.firstvarindex + vi);
  yaml << YAML::EndSeq;
  yaml << YAML::Key << "do_checkpoint" << YAML::Value
       << commongroupdata.do_checkpoint;
  yaml << YAML::Key << "do_restrict" << YAML::Value
       << commongroupdata.do_restrict;
  yaml << YAML::Key << "valid" << YAML::Value << commongroupdata.valid;
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &
operator<<(YAML::Emitter &yaml,
           const GHExt::GlobalData::ArrayGroupData &arraygroupdata) {
  yaml << YAML::LocalTag("arraygroupdata-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "commongroupdata" << YAML::Value
       << (GHExt::CommonGroupData)arraygroupdata;
  yaml << YAML::Key << "data" << YAML::Value << arraygroupdata.data;
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml,
                          const GHExt::GlobalData &globaldata) {
  yaml << YAML::LocalTag("globaldata-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "arraygroupdata" << YAML::Value << YAML::BeginSeq;
  for (const auto &arraygroupdata : globaldata.arraygroupdata)
    if (arraygroupdata)
      yaml << *arraygroupdata;
  yaml << YAML::EndSeq;
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml,
                          const GHExt::LevelData::GroupData &groupdata) {
  yaml << YAML::LocalTag("groupdata-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "commongroupdata" << YAML::Value
       << (GHExt::CommonGroupData)groupdata;
  yaml << YAML::Key << "level" << YAML::Value << groupdata.level;
  yaml << YAML::Key << "indextype" << YAML::Value << YAML::Flow
       << seq(groupdata.indextype);
  yaml << YAML::Key << "nghostzones" << YAML::Value << YAML::Flow
       << seq(groupdata.nghostzones);
  yaml << YAML::Key << "parities" << YAML::Value << YAML::Flow
       << seqs(groupdata.parities);
  yaml << YAML::Key << "fluxes" << YAML::Value << YAML::Flow << YAML::BeginSeq;
  for (const int flux : groupdata.fluxes)
    yaml << (flux >= 0 ? CCTK_FullGroupName(flux) : "");
  yaml << YAML::EndSeq;
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml,
                          const GHExt::LevelData &leveldata) {
  yaml << YAML::LocalTag("leveldata-1.0.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "level" << YAML::Value << leveldata.level;
  yaml << YAML::Key << "is_subcycling_level" << YAML::Value
       << leveldata.is_subcycling_level;
  yaml << YAML::Key << "iteration" << YAML::Value << leveldata.iteration;
  yaml << YAML::Key << "delta_iteration" << YAML::Value
       << leveldata.delta_iteration;
  yaml << YAML::Key << "fab" << YAML::Value << *leveldata.fab;
  yaml << YAML::Key << "groupdata" << YAML::Value << YAML::BeginSeq;
  for (const auto &groupdata : leveldata.groupdata)
    if (groupdata)
      yaml << *groupdata;
  yaml << YAML::EndSeq;
  yaml << YAML::EndMap;
  return yaml;
}

YAML::Emitter &operator<<(YAML::Emitter &yaml, const GHExt &ghext) {
  yaml << YAML::LocalTag("ghext-1.1.0");
  yaml << YAML::BeginMap;
  yaml << YAML::Key << "amrcore" << YAML::Value << *ghext.amrcore;
  yaml << YAML::Key << "globaldata" << YAML::Value << ghext.globaldata;
  yaml << YAML::Key << "leveldata" << YAML::Value << ghext.leveldata;
  yaml << YAML::EndMap;
  return yaml;
}

ostream &operator<<(ostream &os, const GHExt &ghext) {
  YAML::Emitter yaml;
  yaml << ghext;
  os << yaml.c_str() << "\n";
  return os;
}

////////////////////////////////////////////////////////////////////////////////

// Start driver
extern "C" int CarpetX_Startup() {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("Startup");

  // Output a startup message
  const vector<string> features{
#ifdef AMREX_USE_MPI
      "MPI",
#else
      "no MPI",
#endif
#ifdef AMREX_USE_OMP
      "OpenMP",
#else
      "no OpenMP",
#endif
#ifdef AMREX_USE_GPU
      "accelerators",
#else
      "no accelerators",
#endif
#ifdef AMREX_USE_ASSERTION
      "debug",
#else
      "optimized",
#endif
  };
  ostringstream buf;
  buf << logo();
  buf << "AMR driver provided by CarpetX,\n"
      << "using AMReX " << amrex::Version() << " (";
  auto sep = "";
  for (const auto &feature : features) {
    buf << sep << feature;
    sep = ", ";
  }
  buf << ")";
  int ierr = CCTK_RegisterBanner(buf.str().c_str());
  assert(!ierr);

  // Copy parameters
  Loop::CarpetX_poison_undefined_values = poison_undefined_values;

  // Register a GH extension
  ghext_handle = CCTK_RegisterGHExtension("CarpetX");
  assert(ghext_handle >= 0);
  int iret = CCTK_RegisterGHExtensionSetupGH(ghext_handle, SetupGH);
  assert(iret);
  iret = CCTK_RegisterGHExtensionInitGH(ghext_handle, InitGH);
  assert(iret);
  iret = CCTK_RegisterGHExtensionScheduleTraverseGH(ghext_handle,
                                                    ScheduleTraverseGH);
  assert(iret);

  CCTK_OverloadInitialise(Initialise);
  CCTK_OverloadEvolve(Evolve);
  CCTK_OverloadShutdown(Shutdown);
  CCTK_OverloadOutputGH(OutputGH);

  CCTK_OverloadMyProc(MyProc);
  CCTK_OverloadnProcs(nProcs);
  CCTK_OverloadExit(Exit);
  CCTK_OverloadAbort(Abort);
  CCTK_OverloadBarrier(Barrier);

  CCTK_OverloadEnableGroupStorage(EnableGroupStorage);
  CCTK_OverloadDisableGroupStorage(DisableGroupStorage);
  CCTK_OverloadGroupStorageIncrease(GroupStorageIncrease);
  CCTK_OverloadGroupStorageDecrease(GroupStorageDecrease);

  CCTK_OverloadSyncGroupsByDirI(SyncGroupsByDirI);

  CCTK_OverloadInterpGridArrays(CarpetX_InterpGridArrays);
  CCTK_OverloadGroupDynamicData(GroupDynamicData);
  return 0;
}

// Set up GH extension
void *SetupGH(tFleshConfig *fc, int convLevel, cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("SetupGH");

  assert(fc);
  assert(convLevel == 0);
  assert(cctkGH);

  // Initialize AMReX
  amrex::ParmParse pp;
  // Don't catch Unix signals. If signals are caught, we don't get
  // core files.
  pp.add("amrex.signal_handling", 0);
  // Throw exceptions for failing AMReX assertions. With exceptions,
  // we get core files.
  pp.add("amrex.throw_exception", 1);
  // Set tile size
  pp.addarr("fabarray.mfiter_tile_size",
            vector<int>{max_tile_size_x, max_tile_size_y, max_tile_size_z});
  pamrex = amrex::Initialize(MPI_COMM_WORLD);

  // Create grid structure
  ghext = make_unique<GHExt>();

  return ghext.get();
}

// Initialize GH extension
int InitGH(cGH *restrict cctkGH) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("InitGH");

  assert(cctkGH);

  // Domain
  const amrex::RealBox domain({xmin, ymin, zmin}, {xmax, ymax, zmax});

  // Number of coarse grid cells
  const amrex::Vector<int> ncells{ncells_x, ncells_y, ncells_z};

  const int coord = -1; // undefined?

  // Refinement ratios
  const amrex::Vector<amrex::IntVect> reffacts{}; // empty

  // Periodicity
  const amrex::Array<int, dim> is_periodic{
      periodic || periodic_x, periodic || periodic_y, periodic || periodic_z};

  // Set blocking factors via parameter table since AmrMesh needs to
  // know them when its constructor is running, but there are no
  // constructor arguments for them
  amrex::ParmParse pp;
  pp.add("amr.blocking_factor_x", blocking_factor_x);
  pp.add("amr.blocking_factor_y", blocking_factor_y);
  pp.add("amr.blocking_factor_z", blocking_factor_z);
  pp.add("amr.max_grid_size_x", max_grid_size_x);
  pp.add("amr.max_grid_size_y", max_grid_size_y);
  pp.add("amr.max_grid_size_z", max_grid_size_z);
  pp.add("amr.grid_eff", grid_efficiency);

  ghext->amrcore = make_unique<CactusAmrCore>(
      domain, max_num_levels - 1, ncells, coord, reffacts, is_periodic);

  if (verbose) {
#pragma omp critical
    {
      const int maxnumlevels = ghext->amrcore->maxLevel() + 1;
      for (int level = 0; level < maxnumlevels; ++level) {
        CCTK_VINFO("amrex::Geometry level %d:", level);
        cout << ghext->amrcore->Geom(level) << "\n";
      }
    }
  }

  // #pragma omp critical
  // {
  // CCTK_VINFO("amrex::Boxamrex::Array level %d:", level);
  // cout << ghext->amrcore->boxamrex::Array(level) << "\n";
  // CCTK_VINFO("DistributionMap level %d:", level);
  // cout << ghext->amrcore->DistributionMap(level) << "\n";
  // }

  return 0; // unused
}

// Traverse schedule
int ScheduleTraverseGH(cGH *restrict cctkGH, const char *where) {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("ScheduleTraverseGH iteration %d %s", cctkGH->cctk_iteration,
               where);

  int ierr = CCTK_ScheduleTraverse(where, cctkGH, CallFunction);
  //   if (ierr == 2)
  // #pragma omp critical
  //       CCTK_VINFO("Schedule item \"%s\" not found", where);
  assert(ierr == 0 || ierr == 2);

  return 0; // unused
}

// Shut down driver
void ShutdownADIOS2();
void ShutdownOpenPMD();
extern "C" int CarpetX_Shutdown() {
  DECLARE_CCTK_PARAMETERS;

  if (verbose)
#pragma omp critical
    CCTK_VINFO("Shutdown");

  // Shut down ADIOS2
  ShutdownADIOS2();
  ShutdownOpenPMD();

  // Deallocate memory pools
  mempools.reset();

  if (false) {
    // Should we really do this? Cactus's extension handling mechanism
    // becomes inconsistent once extensions have been unregistered.
    int iret = CCTK_UnregisterGHExtension("CarpetX");
    assert(iret == 0);
  }

  // Deallocate grid hierarchy
  ghext = nullptr;

  // Finalize AMReX
  amrex::Finalize(pamrex);
  pamrex = nullptr;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////

int MyProc(const cGH *restrict cctkGH) {
  if (pamrex) {
    return amrex::ParallelDescriptor::MyProc();
  } else {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return rank;
  }
}

int nProcs(const cGH *restrict cctkGH) {
  if (pamrex) {
    return amrex::ParallelDescriptor::NProcs();
  } else {
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    return size;
  }
}

int Exit(cGH *cctkGH, int retval) {
  if (pamrex) {
    amrex::ParallelDescriptor::Abort(retval);
  } else {
    MPI_Abort(MPI_COMM_WORLD, retval);
  }
  return 0; // unreachable
}

int Abort(cGH *cctkGH, int retval) {
  if (pamrex) {
    amrex::ParallelDescriptor::Abort(retval);
  } else {
    MPI_Abort(MPI_COMM_WORLD, retval);
  }
  return 0; // unreachable
}

int Barrier(const cGH *restrict cctkGH) {
  assert(pamrex);
  amrex::ParallelDescriptor::Barrier();
  return 0;
}

void CarpetX_CallScheduleGroup(void *cctkGH_, const char *groupname) {
  cGH *cctkGH = static_cast<cGH *>(cctkGH_);
  static Timer timer("CallScheduleGroup");
  Interval interval(timer);
  int ierr = CCTK_ScheduleTraverse(groupname, cctkGH, CallFunction);
  assert(!ierr);
}

CCTK_INT
CarpetX_GetDomainSpecification(const CCTK_INT size,
                               CCTK_REAL *restrict const physical_min,
                               CCTK_REAL *restrict const physical_max,
                               CCTK_REAL *restrict const interior_min,
                               CCTK_REAL *restrict const interior_max,
                               CCTK_REAL *restrict const exterior_min,
                               CCTK_REAL *restrict const exterior_max,
                               CCTK_REAL *restrict const spacing) {
  assert(size == dim);
  assert(!empty(ghext->leveldata));
  const auto &fab = *ghext->leveldata.at(0).fab;
  const auto &geom = ghext->amrcore->Geom(0);
  for (int d = 0; d < dim; ++d) {
    // the domain itself (physical faces)
    physical_min[d] = geom.ProbDomain().lo(d);
    physical_max[d] = geom.ProbDomain().hi(d);
    // domain without boundary points (last interior cell centre)
    interior_min[d] = geom.ProbDomain().lo(d) + 0.5 * geom.CellSize(d);
    interior_max[d] = geom.ProbDomain().hi(d) - 0.5 * geom.CellSize(d);
    // domain including boundary/ghost points (last boundary/ghost cell centre)
    exterior_min[d] =
        geom.ProbDomain().lo(d) - (fab.nGrow(d) - 0.5) * geom.CellSize(d);
    exterior_max[d] =
        geom.ProbDomain().hi(d) + (fab.nGrow(d) - 0.5) * geom.CellSize(d);
    spacing[d] = geom.CellSize(d);
  }
  return 0;
}

CCTK_INT CarpetX_GetBoundarySizesAndTypes(const void *const cctkGH_,
                                          const CCTK_INT size,
                                          CCTK_INT *restrict const bndsize,
                                          CCTK_INT *restrict const is_ghostbnd,
                                          CCTK_INT *restrict const is_symbnd,
                                          CCTK_INT *restrict const is_physbnd) {
  DECLARE_CCTK_PARAMETERS;

  const cGH *restrict const cctkGH = static_cast<const cGH *>(cctkGH_);
  assert(size == 2 * dim);

  const array<array<bool, 3>, 2> is_periodic{{
      {{
          periodic || periodic_x,
          periodic || periodic_y,
          periodic || periodic_z,
      }},
      {{
          periodic || periodic_x,
          periodic || periodic_y,
          periodic || periodic_z,
      }},
  }};
  const array<array<bool, 3>, 2> is_reflect{{
      {{
          !!reflection_x,
          !!reflection_y,
          !!reflection_z,
      }},
      {{
          !!reflection_upper_x,
          !!reflection_upper_y,
          !!reflection_upper_z,
      }},
  }};

  for (int d = 0; d < dim; ++d) {
    for (int f = 0; f < 2; ++f) {
      bndsize[2 * d + f] = cctkGH->cctk_nghostzones[d];
      is_ghostbnd[2 * d + f] = cctkGH->cctk_bbox[2 * d + f];
      is_symbnd[2 * d + f] =
          !is_ghostbnd[2 * d + f] && (is_periodic[f][d] || is_reflect[f][d]);
      is_physbnd[2 * d + f] = !is_ghostbnd[2 * d + f] && !is_symbnd[2 * d + f];
    }
  }
  return 0;
}

int GroupDynamicData(const cGH *cctkGH, int gi, cGroupDynamicData *data) {
  // Return values:
  //  0 for success
  // -1 if given pointer to data structure is NULL
  // -3 if given GH pointer is invalid
  // (-77 if group has zero variables)
  // -78 if group does not exist
  if (not cctkGH)
    return -3;
  if (not(gi >= 0 and gi < CCTK_NumGroups()))
    return -78;
  if (not data)
    return -1;
  cGroup group;
  int ierr = CCTK_GroupData(gi, &group);
  assert(!ierr);
  if (group.grouptype == CCTK_GF) {
    data->dim = group.dim;
    data->lsh = cctkGH->cctk_lsh;
    data->ash = cctkGH->cctk_ash;
    data->gsh = cctkGH->cctk_gsh;
    data->lbnd = cctkGH->cctk_lbnd;
    data->ubnd = cctkGH->cctk_ubnd;
    data->tile_min = cctkGH->cctk_tile_min;
    data->tile_max = cctkGH->cctk_tile_max;
    data->bbox = cctkGH->cctk_bbox;
    data->nghostzones = cctkGH->cctk_nghostzones;
    data->activetimelevels = CCTK_ActiveTimeLevelsGI(cctkGH, gi);
  } else if (group.grouptype == CCTK_SCALAR or group.grouptype == CCTK_ARRAY) {
    GHExt::GlobalData &globaldata = ghext->globaldata;
    GHExt::GlobalData::ArrayGroupData &arraygroupdata =
        *globaldata.arraygroupdata.at(gi);
    data->dim = arraygroupdata.dimension;
    data->lsh = arraygroupdata.lsh;
    data->ash = arraygroupdata.ash;
    data->gsh = arraygroupdata.gsh;
    data->lbnd = arraygroupdata.lbnd;
    data->ubnd = arraygroupdata.ubnd;
    data->bbox = arraygroupdata.bbox;
    data->nghostzones = arraygroupdata.nghostzones;
    data->activetimelevels = arraygroupdata.activetimelevels;
  } else {
    CCTK_VERROR("Internal error: unexpected group type %d for group '%s'",
                (int)group.grouptype, CCTK_FullGroupName(gi));
  }
  return 0;
}

} // namespace CarpetX
