#include "io_openpmd.hxx"

#include "driver.hxx"
#include "timer.hxx"

#include <div.hxx>
#include <vect.hxx>

#include <fixmath.hxx>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#ifdef HAVE_CAPABILITY_openPMD_api

#include <nlohmann/json.hpp>
#include <openPMD/openPMD.hpp>

#ifdef _OPENMP
#include <omp.h>
#else
static inline int omp_get_max_threads() { return 1; }
static inline int omp_get_num_threads() { return 1; }
static inline int omp_get_thread_num() { return 0; }
static inline int omp_in_parallel() { return 0; }
#endif

#include <mpi.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstring>
#include <ctime>
#include <ios>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace CarpetX {

constexpr bool io_verbose = true;

constexpr openPMD::IterationEncoding iterationEncoding =
    openPMD::IterationEncoding::fileBased;
// constexpr openPMD::IterationEncoding iterationEncoding =
//     openPMD::IterationEncoding::variableBased;

static constexpr const char suffix[] = "bp"; // ["bp", "json", "h5"]

//  constexpr const char options[]
const nlohmann::json json_options{
    {"adios2",
     {
         {"engine",
          {
              {"type", "BP4"},
              {"parameters",
               {
                   {"BufferGrowthFactor", "2.0"},
                   // {"MaxBufferSize", "500 MB"},
               }},
          }},
         {"dataset",
          {
              {"operators",
               {
                   {
                       {"type", "blosc"},
                       {"parameters",
                        {
                            {"clevel", "9"},
                            {"doshuffle", "BLOSC_BITSHUFFLE"},
                        }},
                   },
               }},
          }},
     }},
};
const std::string options = json_options.dump();

constexpr bool input_ghosts = false;
constexpr bool output_ghosts = false;

////////////////////////////////////////////////////////////////////////////////

struct Const {
  // From: CODATA Internationally recommended 2018 values of the
  // Fundamental Physical Constants
  static constexpr CCTK_REAL c = 299792458;         // m s-1
  static constexpr CCTK_REAL G = 6.67430e-11;       // m3 kg-1 s-2
  static constexpr CCTK_REAL M_solar = 1.98847e+30; // kg
};

struct Unit {
  // We use c = G = 1, and M_solar as mass unit.
  static constexpr CCTK_REAL velocity = Const::c;   // m s-1
  static constexpr CCTK_REAL mass = Const::M_solar; // kg
  static constexpr CCTK_REAL length = Const::G * mass / pow(Const::c, 2); // m
  static constexpr CCTK_REAL time = length / velocity;                    // s
};

////////////////////////////////////////////////////////////////////////////////

namespace {

template <typename T, std::size_t N>
constexpr std::array<T, N> reversed(const std::array<T, N> &arr) {
  std::array<T, N> res;
  for (std::size_t n = 0; n < N; ++n)
    res[n] = arr[N - 1 - n];
  return res;
}

template <typename T, int D>
constexpr Arith::vect<T, D> reversed(const Arith::vect<T, D> &arr) {
  Arith::vect<T, D> res;
  for (std::size_t d = 0; d < D; ++d)
    res[d] = arr[D - 1 - d];
  return res;
}

template <typename T> std::vector<T> reversed(const std::vector<T> &vec) {
  std::vector<T> res(vec.size());
  std::reverse_copy(vec.begin(), vec.end(), res.begin());
  return res;
}
template <typename T> std::vector<T> reversed(std::vector<T> &&vec) {
  std::vector<T> res(vec);
  std::reverse(res.begin(), res.end());
  return res;
}

template <typename T = std::uint64_t, typename C>
std::vector<T> to_vector(const C &source) {
  std::vector<T> res(source.size());
  for (std::size_t n = 0; n < res.size(); ++n)
    res[n] = source[n];
  return res;
}

template <typename T> std::string vector_string(const std::vector<T> &vec) {
  std::ostringstream buf;
  buf << "[";
  for (std::size_t n = 0; n < vec.size(); ++n) {
    if (n != 0)
      buf << ",";
    buf << vec[n];
  }
  buf << "]";
  return buf.str();
}

} // namespace

////////////////////////////////////////////////////////////////////////////////

std::string Geometry_string(const openPMD::Mesh::Geometry &geometry) {
  std::ostringstream buf;
  buf << geometry;
  return buf.str();
}

std::string UnitDimension_string(const std::array<double, 7> &unitDimension) {
  std::ostringstream buf;
  buf << "L:" << unitDimension[uint8_t(openPMD::UnitDimension::L)] << " "
      << "M:" << unitDimension[uint8_t(openPMD::UnitDimension::M)] << " "
      << "T:" << unitDimension[uint8_t(openPMD::UnitDimension::T)] << " "
      << "I:" << unitDimension[uint8_t(openPMD::UnitDimension::I)] << " "
      << "θ:" << unitDimension[uint8_t(openPMD::UnitDimension::theta)] << " "
      << "N:" << unitDimension[uint8_t(openPMD::UnitDimension::N)] << " "
      << "J:" << unitDimension[uint8_t(openPMD::UnitDimension::J)];
  return buf.str();
}

////////////////////////////////////////////////////////////////////////////////

struct carpetx_openpmd_t {
  carpetx_openpmd_t() = default;

  carpetx_openpmd_t(const carpetx_openpmd_t &) = delete;
  carpetx_openpmd_t(carpetx_openpmd_t &&) = default;
  carpetx_openpmd_t &operator=(const carpetx_openpmd_t &) = delete;
  carpetx_openpmd_t &operator=(carpetx_openpmd_t &&) = default;

  static std::optional<carpetx_openpmd_t> self;

  ////////////////////////////////////////////////////////////////////////////////

  template <typename T, std::size_t D> struct box_t {
    Arith::vect<T, D> lo, hi;
    constexpr Arith::vect<T, D> shape() const {
      Arith::vect<T, D> sh;
      for (std::size_t d = 0; d < D; ++d)
        sh[d] = hi[d] < lo[d] ? T{0} : hi[d] - lo[d];
      return sh;
    }
    constexpr T size() const {
      const auto sh = shape();
      T sz{1};
      for (std::size_t d = 0; d < D; ++d)
        sz *= sh[d];
      return sz;
    }
    constexpr Arith::vect<T, D> stride() const {
      const Arith::vect<T, D> sh = shape();
      Arith::vect<T, D> str;
      T np{1};
      for (std::size_t d = 0; d < D; ++d)
        str[d] = (np *= sh[d]);
      return str;
    }
    constexpr T linear(const Arith::vect<T, D> &index) const {
      const Arith::vect<T, D> str = stride();
      T lin{0};
      if (D > 0) {
        lin = index[0];
        for (std::size_t d = 1; d < D; ++d)
          lin += index[d] * str[d - 1];
      }
      return lin;
    }
  };

  template <typename T, typename I, std::size_t D> struct level_t {
    box_t<T, D> rdomain;
    Arith::vect<bool, D> is_cell_centred;
    box_t<I, D> idomain;
    constexpr Arith::vect<T, D> rcoord(const Arith::vect<I, D> &icoord) const {
      Arith::vect<T, D> r;
      for (std::size_t d = 0; d < D; ++d) {
        const T rlo = rdomain.lo[d];
        const T rhi = rdomain.hi[d];
        const I ilo2 = 2 * idomain.lo[d] - is_cell_centred;
        const I ihi2 = 2 * idomain.hi[d] + is_cell_centred;
        const I i2 = 2 * icoord[d];
        // The expression below has been carefully constructed to
        // avoid round-off errors at i=ilo and i=ihi. Do not rearrange
        // the terms without preserving this property.
        r[d] = (T(ihi2 - i2) / T(ihi2 - ilo2)) * rlo[d] +
               (T(i2 - ilo2) / T(ihi2 - ilo2)) * rhi[d];
      }
      return r;
    }
    std::vector<box_t<I, D> > grids;
    Arith::vect<std::vector<I>, 2> offsets_sizes() const {
      std::vector<I> offsets(grids.size() + 1), sizes(grids.size());
      I offset{0};
      for (std::size_t n = 0; n < grids.size(); ++n) {
        const T size = grids[n].size();
        offsets[n] = offset;
        sizes[n] = size;
        offset += size;
      }
      offsets[grids.size()] = offset;
      return {offsets, sizes};
    }
  };

  template <typename T, typename I, std::size_t D> struct grid_structure_t {
    box_t<T, D> rdomain;
    std::vector<level_t<T, I, D> > levels;
  };

  ////////////////////////////////////////////////////////////////////////////////

  // Allowed characters are only [A-Za-z_]
  static std::string make_meshname(const int gi, const int level) {
    std::string groupname = CCTK_FullGroupName(gi);
    groupname = std::regex_replace(groupname, std::regex("::"), "_");
    for (auto &ch : groupname)
      ch = std::tolower(ch);
    std::ostringstream buf;
    buf << groupname;
    buf << "_rl" << setw(2) << setfill('0') << level;
    return buf.str();
  }

#if 0
  static std::tuple<int, int> interpret_meshname(const std::string &meshname) {
    std::smatch match;
    const bool matched =
        std::regex_match(meshname, match, std::regex("(\\w+)_rl0*(\\d+)"));
    if (!matched)
      CCTK_VERROR("Cannot parse mesh name %s", meshname.c_str());
    const std::string groupname = match[1].str();
    const int gi = CCTK_GroupIndex(groupname.c_str());
    if (gi < 0)
      CCTK_VERROR("Unknown group name %s", groupname.c_str());
    const std::string levelstr = match[2].str();
    const int level = std::stoi(levelstr);
    return {gi, level};
  }
#endif

  // Allowed characters are only [A-Za-z_]
  static std::string make_componentname(const int gi, const int vi) {
    const int v0 = CCTK_FirstVarIndexI(gi);
    std::string varname = CCTK_FullVarName(v0 + vi);
    varname = std::regex_replace(varname, std::regex("::"), "_");
    for (auto &ch : varname)
      ch = std::tolower(ch);
    return varname;
  }

  ////////////////////////////////////////////////////////////////////////////////

  std::optional<std::string> filename;
  std::optional<openPMD::Series> series;
  std::optional<openPMD::WriteIterations> write_iters;

  void InputOpenPMD(const cGH *const cctkGH,
                    const std::vector<bool> &input_group,
                    const std::string &input_dir,
                    const std::string &input_file);

  void OutputOpenPMD(const cGH *const cctkGH,
                     const std::vector<bool> &output_group,
                     const std::string &output_dir,
                     const std::string &output_file);
};

////////////////////////////////////////////////////////////////////////////////

std::optional<carpetx_openpmd_t> carpetx_openpmd_t::self;

void InputOpenPMD(const cGH *cctkGH, const std::vector<bool> &input_group,
                  const std::string &input_dir, const std::string &input_file) {
  if (!carpetx_openpmd_t::self)
    carpetx_openpmd_t::self = std::make_optional<carpetx_openpmd_t>();
  carpetx_openpmd_t::self->InputOpenPMD(cctkGH, input_group, input_dir,
                                        input_file);
}

void OutputOpenPMD(const cGH *const cctkGH,
                   const std::vector<bool> &output_group,
                   const std::string &output_dir,
                   const std::string &output_file) {
  if (!carpetx_openpmd_t::self)
    carpetx_openpmd_t::self = std::make_optional<carpetx_openpmd_t>();
  carpetx_openpmd_t::self->OutputOpenPMD(cctkGH, output_group, output_dir,
                                         output_file);
}

void ShutdownOpenPMD() { carpetx_openpmd_t::self.reset(); }

////////////////////////////////////////////////////////////////////////////////

void carpetx_openpmd_t::InputOpenPMD(const cGH *const cctkGH,
                                     const std::vector<bool> &input_group,
                                     const std::string &input_dir,
                                     const std::string &input_file) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const int lapse_gi = CCTK_GroupIndex("ADMBase::lapse");
  assert(lapse_gi >= 0);

  // Set up timers
  static Timer timer("InputOpenPMD");
  Interval interval(timer);

  if (std::count(input_group.begin(), input_group.end(), true) == 0)
    return;

  const bool is_root = CCTK_MyProc(nullptr) == 0;
  if (is_root) {
    CCTK_VINFO("openPMD input for groups:");
    for (int gi = 0; gi < CCTK_NumGroups(); ++gi)
      if (input_group.at(gi))
        CCTK_VINFO("  %s", CCTK_FullGroupName(gi));
  }

  if (io_verbose)
    CCTK_VINFO("InputOpenPMD...");

  if (!series) {

    if (io_verbose)
      CCTK_VINFO("Creating openPMD object...");
    std::ostringstream buf;
    switch (iterationEncoding) {
    case openPMD::IterationEncoding::fileBased:
      buf << input_dir << "/" << input_file << ".it%08T." << suffix;
      break;
    case openPMD::IterationEncoding::variableBased:
      buf << input_dir << "/" << input_file << "." << suffix;
      break;
    default:
      abort();
    }
    filename = std::make_optional<std::string>(buf.str());
    series = std::make_optional<openPMD::Series>(
        *filename, openPMD::Access::READ_ONLY, MPI_COMM_WORLD, options);
  }
  assert(filename);
  assert(series);

  openPMD::ReadIterations readIters = series->readIterations();
  assert(readIters.begin() != readIters.end());
  openPMD::IndexedIteration iter = *readIters.begin();
  const uint64_t iterIndex = iter.iterationIndex;
  CCTK_VINFO("  iteration: %d", int(iterIndex));

  const CCTK_REAL time = iter.time<CCTK_REAL>();
  const CCTK_REAL dt = iter.dt<CCTK_REAL>();
  const double timeUnitSI = iter.timeUnitSI();
  CCTK_VINFO("  time: %f", double(time));
  CCTK_VINFO("  dt: %f", double(dt));
  CCTK_VINFO("  time unit SI: %f", timeUnitSI);

  openPMD::Container<openPMD::Mesh> &meshes = iter.meshes;
  CCTK_VINFO("  found %d meshes", int(meshes.size()));

#if 1
  for (auto mesh_iter = meshes.begin(); mesh_iter != meshes.end();
       ++mesh_iter) {
    const std::string &mesh_name = mesh_iter->first;
    openPMD::Mesh &mesh = mesh_iter->second;
    CCTK_VINFO("    mesh: %s", mesh_name.c_str());

    const openPMD::Mesh::Geometry geometry = mesh.geometry();
    const std::vector<std::string> axisLabels = mesh.axisLabels();
    const std::vector<CCTK_REAL> gridSpacing = mesh.gridSpacing<CCTK_REAL>();
    const std::vector<double> gridGlobalOffset = mesh.gridGlobalOffset();
    const double gridUnitSI = mesh.gridUnitSI();
    const std::array<double, 7> unitDimension = mesh.unitDimension();
    const CCTK_REAL timeOffset = mesh.timeOffset<CCTK_REAL>();
    CCTK_VINFO("      geometry: %s", Geometry_string(geometry).c_str());
    CCTK_VINFO("      axis labels: %s",
               vector_string(reversed(axisLabels)).c_str());
    CCTK_VINFO("      grid spacing: %s",
               vector_string(reversed(gridSpacing)).c_str());
    CCTK_VINFO("      grid global offset: %s",
               vector_string(reversed(gridGlobalOffset)).c_str());
    CCTK_VINFO("      grid unit SI: %f", gridUnitSI);
    CCTK_VINFO("      unit dimension: %s",
               UnitDimension_string(unitDimension).c_str());
    CCTK_VINFO("      time offset: %f", double(timeOffset));

    CCTK_VINFO("      found %d components", int(mesh.size()));
    openPMD::Extent extent;
    for (auto record_component_iter = mesh.begin();
         record_component_iter != mesh.end(); ++record_component_iter) {
      const std::string &record_component_name = record_component_iter->first;
      openPMD::MeshRecordComponent &record_component =
          record_component_iter->second;
      CCTK_VINFO("        component: %s", record_component_name.c_str());

      const std::vector<CCTK_REAL> position =
          record_component.position<CCTK_REAL>();
      CCTK_VINFO("          position: %s",
                 vector_string(reversed(position)).c_str());

      const int ndims = record_component.getDimensionality();
      if (extent.empty())
        extent = record_component.getExtent();
      else
        assert(extent == record_component.getExtent());
      CCTK_VINFO("          ndims: %d", ndims);
      CCTK_VINFO("          extent: %s",
                 vector_string(reversed(extent)).c_str());

      const std::vector<openPMD::WrittenChunkInfo> chunks =
          record_component.availableChunks();
      CCTK_VINFO("          found %d chunks", int(chunks.size()));
      if (mesh_iter == meshes.begin() &&
          record_component_iter == mesh.begin()) {
        for (std::size_t n = 0; n < chunks.size(); ++n) {
          const openPMD::WrittenChunkInfo &chunk = chunks.at(n);
          CCTK_VINFO("            chunk: %d   start: %s   count: %s", int(n),
                     vector_string(reversed(chunk.offset)).c_str(),
                     vector_string(reversed(chunk.extent)).c_str());
        }
      }

    } // for record_component
  }   // for mesh
#endif

  // Loop over levels
  for (const auto &leveldata : ghext->leveldata) {
    if (io_verbose)
      CCTK_VINFO("Reading level %d", leveldata.level);

    // Determine grid structure

    const int *const nghosts = cctkGH->cctk_nghostzones;
    const amrex::Geometry &geom = ghext->amrcore->Geom(leveldata.level);
    const double *const xlo = geom.ProbLo();
    const double *const xhi = geom.ProbHi();
    const double *const dx = geom.CellSize();
    const box_t<CCTK_REAL, 3> rdomain{
      lo : {xlo[0] - nghosts[0] * dx[0], xlo[1] - nghosts[1] * dx[1],
            xlo[2] - nghosts[2] * dx[2]},
      hi : {xhi[0] + nghosts[0] * dx[0], xhi[1] + nghosts[1] * dx[1],
            xhi[2] + nghosts[2] * dx[2]}
    };
    const amrex::Box &dom = geom.Domain();
    const amrex::IntVect &ilo = dom.smallEnd();
    const amrex::IntVect &ihi = dom.bigEnd();
    // The domain is always vertex centred. The tensor components are
    // then staggered if necessary.
    const box_t<int, 3> idomain{
      lo : {ilo[0] - nghosts[0], ilo[1] - nghosts[1], ilo[2] - nghosts[2]},
      hi : {ihi[0] + nghosts[0] + 1 + 1, ihi[1] + nghosts[1] + 1 + 1,
            ihi[2] + nghosts[2] + 1 + 1}
    };
    if (io_verbose) {
      CCTK_VINFO("Level: %d", leveldata.level);
      CCTK_VINFO("  xmin: [%f,%f,%f]", double(rdomain.lo[0]),
                 double(rdomain.lo[1]), double(rdomain.lo[2]));
      CCTK_VINFO("  xmax: [%f,%f,%f]", double(rdomain.hi[0]),
                 double(rdomain.hi[1]), double(rdomain.hi[2]));
      CCTK_VINFO("  imin: [%d,%d,%d]", int(idomain.lo[0]), int(idomain.lo[1]),
                 int(idomain.lo[2]));
      CCTK_VINFO("  imax: [%d,%d,%d]", int(idomain.hi[0]), int(idomain.hi[1]),
                 int(idomain.hi[2]));
    }

    const int numgroups = CCTK_NumGroups();
    for (int gi = 0; gi < numgroups; ++gi) {
      if (input_group.at(gi)) {
        if (io_verbose)
          CCTK_VINFO("Reading group %d %s...", gi, CCTK_FullGroupName(gi));

        // Check group properties

        cGroup cgroup;
        const int ierr = CCTK_GroupData(gi, &cgroup);
        assert(!ierr);
        assert(cgroup.grouptype == CCTK_GF);
        assert(cgroup.vartype == CCTK_VARIABLE_REAL);
        assert(cgroup.dim == 3);
        // cGroupDynamicData cgroupdynamicdata;
        // ierr = CCTK_GroupDynamicData(cctkGH, gi, &cgroupdynamicdata);
        // assert(!ierr);
        // TODO: Check whether group has storage
        // TODO: Check whether data are valid

        auto &groupdata = *leveldata.groupdata.at(gi);
        // const int firstvarindex = groupdata.firstvarindex;
        const int numvars = groupdata.numvars;
        const int tl = 0;
        amrex::MultiFab &mfab = *groupdata.mfab[tl];
        const amrex::IndexType &indextype = mfab.ixType();
        const Arith::vect<bool, 3> is_cell_centred{indextype.cellCentered(0),
                                                   indextype.cellCentered(1),
                                                   indextype.cellCentered(2)};

        const int num_local_components = mfab.local_size();

        // Read mesh

        const std::string meshname = make_meshname(gi, leveldata.level);
        if (io_verbose)
          CCTK_VINFO("Reading mesh %s...", meshname.c_str());
        assert(iter.meshes.count(meshname));
        const openPMD::Mesh &mesh = iter.meshes.at(meshname);

        // Define tensor components

        // TODO: Set component names according to the tensor type
        std::vector<openPMD::MeshRecordComponent> record_components;
        record_components.reserve(numvars);
        openPMD::Extent extent;
        for (int vi = 0; vi < numvars; ++vi) {
          const std::string componentname = make_componentname(gi, vi);
          assert(mesh.count(componentname));
          record_components.push_back(mesh.at(componentname));
          const openPMD::MeshRecordComponent &record_component =
              record_components.back();
          if (vi == 0)
            extent = record_component.getExtent();
          else
            assert(extent == record_component.getExtent());
        }
        assert(int(record_components.size()) == numvars);

        // Read data

        if (io_verbose)
          CCTK_VINFO("Reading %d variables with %d components...", numvars,
                     num_local_components);

        // Loop over components (AMReX boxes)
        for (int local_component = 0; local_component < num_local_components;
             ++local_component) {
          const int component = mfab.IndexArray().at(local_component);

          const amrex::Box &fabbox =
              mfab.fabbox(component); // exterior (with ghosts)
          const box_t<int, 3> extbox{
            lo : {fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)},
            hi : {fabbox.bigEnd(0) + 1, fabbox.bigEnd(1) + 1,
                  fabbox.bigEnd(2) + 1}
          };
          const amrex::Box &validbox =
              mfab.box(component); // interior (without ghosts)
          const box_t<int, 3> intbox{
            lo : {validbox.smallEnd(0), validbox.smallEnd(1),
                  validbox.smallEnd(2)},
            hi : {validbox.bigEnd(0) + 1, validbox.bigEnd(1) + 1,
                  validbox.bigEnd(2) + 1}
          };
          const box_t<int, 3> &box = output_ghosts ? extbox : intbox;

          const openPMD::Offset start =
              to_vector(reversed(box.lo - idomain.lo));
          const openPMD::Extent count = to_vector(reversed(box.shape()));
          const int np = box.size();
          assert(int(count.at(0) * count.at(1) * count.at(2)) == np);
          for (int d = 0; d < 3; ++d)
            assert(start.at(d) >= 0);
          for (int d = 0; d < 3; ++d)
            assert(start.at(d) + count.at(d) <= extent.at(d));

          amrex::FArrayBox &fab = mfab[component];
          for (int vi = 0; vi < numvars; ++vi) {

            if (input_ghosts) {
              CCTK_REAL *const ptr = fab.dataPtr() + vi * np;
              record_components.at(vi).loadChunk(openPMD::shareRaw(ptr), start,
                                                 count);
            } else {
              const int amrex_size = extbox.size();
              CCTK_REAL *const ptr = fab.dataPtr() + vi * amrex_size;
              const Arith::vect<int, 3> amrex_shape = extbox.shape();
              const Arith::vect<int, 3> amrex_offset = box.lo - extbox.lo;
              constexpr int amrex_di = 1;
              const int amrex_dj = amrex_di * amrex_shape[0];
              const int amrex_dk = amrex_dj * amrex_shape[1];
              // const int amrex_np = amrex_dk * amrex_shape[2];
              CCTK_REAL *const amrex_ptr = ptr + amrex_di * amrex_offset[0] +
                                           amrex_dj * amrex_offset[1] +
                                           amrex_dk * amrex_offset[2];
              const Arith::vect<int, 3> contig_shape = box.shape();
              constexpr int contig_di = 1;
              const int contig_dj = contig_di * contig_shape[0];
              const int contig_dk = contig_dj * contig_shape[1];
              const int contig_np = contig_dk * contig_shape[2];
              assert(contig_np == np);
              CCTK_REAL *const contig_ptr = ptr + extbox.size() - box.size();
              // TODO: optimize memory layout
#if 0
              CCTK_REAL *const contig_ptr =
                  ptr + contig_di * (contig_shape[0] - 1) +
                  contig_dj * (contig_shape[1] - 1) +
                  contig_dk * (contig_shape[2] - 1) + 1 - contig_np;
              assert(&amrex_ptr[amrex_di * (amrex_shape[0] - 1) +
                                amrex_dj * (amrex_shape[1] - 1) +
                                amrex_dk * (amrex_shape[2] - 1)] ==
                     &contig_ptr[contig_di * (contig_shape[0] - 1) +
                                 contig_dj * (contig_shape[1] - 1) +
                                 contig_dk * (contig_shape[2] - 1)]);
#endif
              const int reflevel = leveldata.level;
              const auto expand_box = [=](const CCTK_REAL
                                              *restrict const contig_ptr) {
                if (gi == lapse_gi)
                  CCTK_VINFO("Expanding lapse on level %d component %d",
                             reflevel, local_component);
                for (int k = 0; k < contig_shape[2]; ++k)
                  for (int j = 0; j < contig_shape[1]; ++j)
                    for (int i = 0; i < contig_shape[0]; ++i)
                      amrex_ptr[amrex_di * i + amrex_dj * j + amrex_dk * k] =
                          contig_ptr[contig_di * i + contig_dj * j +
                                     contig_dk * k];
                if (gi == lapse_gi) {
                  CCTK_VINFO("amrex_shape=[%d,%d,%d]", amrex_shape[0],
                             amrex_shape[1], amrex_shape[2]);
                  CCTK_VINFO("contig_shape=[%d,%d,%d]", contig_shape[0],
                             contig_shape[1], contig_shape[2]);
                  for (int k = 0; k < contig_shape[2]; ++k)
                    for (int j = 0; j < contig_shape[1]; ++j)
                      for (int i = 0; i < contig_shape[0]; ++i)
                        if (amrex_ptr[amrex_di * i + amrex_dj * j +
                                      amrex_dk * k] == 0 ||
                            !isfinite(amrex_ptr[amrex_di * i + amrex_dj * j +
                                                amrex_dk * k]))
                          CCTK_VINFO(
                              "alpha[rl=%d; %d,%d,%d]=%f", reflevel, i, j, k,
                              double(amrex_ptr[amrex_di * i + amrex_dj * j +
                                               amrex_dk * k]));
                }
              };
              if (gi == lapse_gi)
                CCTK_VINFO("Scheduling reading lapse on level %d component %d",
                           reflevel, local_component);
              if (gi == lapse_gi)
                for (int k = 0; k < contig_shape[2]; ++k)
                  for (int j = 0; j < contig_shape[1]; ++j)
                    for (int i = 0; i < contig_shape[0]; ++i)
                      ptr[i * amrex_di + j * amrex_dj + k * amrex_dk] = 0;
              record_components.at(vi).loadChunk(
                  std::shared_ptr<CCTK_REAL>(contig_ptr, expand_box), start,
                  count);
            }

            // Mark read variables as valid
            groupdata.valid.at(tl).at(vi).set(
                input_ghosts ? make_valid_all() : make_valid_int(),
                []() { return "read from openPMD file"; });

          } // for vi
        }   // for local_component
      }
    } // for gi

  } // for leveldata

  if (io_verbose)
    CCTK_VINFO("InputOpenPMD: Performing all reads...");
  series->flush();

  if (io_verbose)
    CCTK_VINFO("InputOpenPMD: Closing iteration...");
  iter.close();

  if (io_verbose)
    CCTK_VINFO("InputOpenPMD: Deallocating objects...");
  series.reset();
  filename.reset();

  if (io_verbose)
    CCTK_VINFO("InputOpenPMD done.");

  if (io_verbose)
    timer.print();
}

////////////////////////////////////////////////////////////////////////////////

void carpetx_openpmd_t::OutputOpenPMD(const cGH *const cctkGH,
                                      const std::vector<bool> &output_group,
                                      const std::string &output_dir,
                                      const std::string &output_file) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  // Set up timers
  static Timer timer("OutputOpenPMD");
  Interval interval(timer);

  if (std::count(output_group.begin(), output_group.end(), true) == 0)
    return;

  if (io_verbose)
    CCTK_VINFO("OutputOpenPMD...");

  if (!series) {

    if (io_verbose)
      CCTK_VINFO("Creating openPMD object...");
    const int mode = 0755;
    static once_flag create_directory;
    call_once(create_directory, [&]() {
      const int ierr = CCTK_CreateDirectory(mode, output_dir.c_str());
      assert(ierr >= 0);
    });
    std::ostringstream buf;
    switch (iterationEncoding) {
    case openPMD::IterationEncoding::fileBased:
      buf << output_dir << "/" << output_file << ".it%08T." << suffix;
      break;
    case openPMD::IterationEncoding::variableBased:
      buf << output_dir << "/" << output_file << "." << suffix;
      break;
    default:
      abort();
    }
    filename = std::make_optional<std::string>(buf.str());
    static bool is_first_output = true;
    const openPMD::Access access =
        iterationEncoding == openPMD::IterationEncoding::fileBased
            ? openPMD::Access::CREATE
        : iterationEncoding == openPMD::IterationEncoding::variableBased
            ? (is_first_output ? openPMD::Access::CREATE
                               : openPMD::Access::READ_WRITE)
            : openPMD::Access::READ_ONLY /*error*/;
    is_first_output = false;
    CCTK_VINFO("  options: %s", options.c_str());
    series = std::make_optional<openPMD::Series>(*filename, access,
                                                 MPI_COMM_WORLD, options);
    series->setIterationEncoding(iterationEncoding);

    series->setAuthor(out_openpmd_author);
    // Software is always "openPMD-api"
    // series->setSoftware("Einstein Toolkit <https://einsteintoolkit.org>");
    // series->setSoftwareVersion("...");
    // Date is set automatically
    // const std::time_t t = std::time(nullptr);
    // char date[100];
    // std::strftime(date, sizeof date, "%Y-%m-%dT%H:%M:%S",
    // std::localtime(&t)); series->setDate(date);
    // series->setSoftwareDependencies("...");
    // series->setMachine("...");

    write_iters =
        std::make_optional<openPMD::WriteIterations>(series->writeIterations());
  }
  assert(filename);
  assert(series);
  assert(write_iters);

  if (io_verbose)
    CCTK_VINFO("Creating iteration %d...", cctk_iteration);
  openPMD::Iteration iter = (*write_iters)[cctk_iteration];
  iter.setTime(cctk_time);
  iter.setDt(cctk_delta_time);
  iter.setTimeUnitSI(Unit::time);

  // Loop over levels
  for (const auto &leveldata : ghext->leveldata) {
    if (io_verbose)
      CCTK_VINFO("Writing level mesh %d...", leveldata.level);

    // Determine grid structure

    const int *const nghosts = cctkGH->cctk_nghostzones;
    const amrex::Geometry &geom = ghext->amrcore->Geom(leveldata.level);
    const double *const xlo = geom.ProbLo();
    const double *const xhi = geom.ProbHi();
    const double *const dx = geom.CellSize();
    const box_t<CCTK_REAL, 3> rdomain{
      lo : {xlo[0] - nghosts[0] * dx[0], xlo[1] - nghosts[1] * dx[1],
            xlo[2] - nghosts[2] * dx[2]},
      hi : {xhi[0] + nghosts[0] * dx[0], xhi[1] + nghosts[1] * dx[1],
            xhi[2] + nghosts[2] * dx[2]}
    };
    const amrex::Box &dom = geom.Domain();
    const amrex::IntVect &ilo = dom.smallEnd();
    const amrex::IntVect &ihi = dom.bigEnd();
    // The domain is always vertex centred. The tensor components are
    // then staggered if necessary.
    const box_t<int, 3> idomain{
      lo : {ilo[0] - nghosts[0], ilo[1] - nghosts[1], ilo[2] - nghosts[2]},
      hi : {ihi[0] + nghosts[0] + 1 + 1, ihi[1] + nghosts[1] + 1 + 1,
            ihi[2] + nghosts[2] + 1 + 1}
    };
    if (io_verbose) {
      CCTK_VINFO("Level: %d", leveldata.level);
      CCTK_VINFO("  xmin: [%f,%f,%f]", double(rdomain.lo[0]),
                 double(rdomain.lo[1]), double(rdomain.lo[2]));
      CCTK_VINFO("  xmax: [%f,%f,%f]", double(rdomain.hi[0]),
                 double(rdomain.hi[1]), double(rdomain.hi[2]));
      CCTK_VINFO("  imin: [%d,%d,%d]", int(idomain.lo[0]), int(idomain.lo[1]),
                 int(idomain.lo[2]));
      CCTK_VINFO("  imax: [%d,%d,%d]", int(idomain.hi[0]), int(idomain.hi[1]),
                 int(idomain.hi[2]));
    }

    // Create dataset

    const openPMD::Datatype datatype = openPMD::determineDatatype<CCTK_REAL>();
    const openPMD::Extent extent = to_vector(reversed(idomain.shape()));
    const openPMD::Dataset dataset(datatype, extent);

    const int numgroups = CCTK_NumGroups();
    for (int gi = 0; gi < numgroups; ++gi) {
      if (output_group.at(gi)) {
        if (io_verbose)
          CCTK_VINFO("Writing group %d %s...", gi, CCTK_FullGroupName(gi));

        // Check group properties

        cGroup cgroup;
        const int ierr = CCTK_GroupData(gi, &cgroup);
        assert(!ierr);
        assert(cgroup.grouptype == CCTK_GF);
        assert(cgroup.vartype == CCTK_VARIABLE_REAL);
        assert(cgroup.dim == 3);
        // cGroupDynamicData cgroupdynamicdata;
        // ierr = CCTK_GroupDynamicData(cctkGH, gi, &cgroupdynamicdata);
        // assert(!ierr);
        // TODO: Check whether group has storage
        // TODO: Check whether data are valid

        const auto &groupdata = *leveldata.groupdata.at(gi);
        // const int firstvarindex = groupdata.firstvarindex;
        const int numvars = groupdata.numvars;
        const int tl = 0;
        const amrex::MultiFab &mfab = *groupdata.mfab[tl];
        const amrex::IndexType &indextype = mfab.ixType();
        const Arith::vect<bool, 3> is_cell_centred{indextype.cellCentered(0),
                                                   indextype.cellCentered(1),
                                                   indextype.cellCentered(2)};

        const int num_local_components = mfab.local_size();

        // Create mesh

        const std::string meshname = make_meshname(gi, leveldata.level);
        if (io_verbose)
          CCTK_VINFO("Defining mesh %s...", meshname.c_str());
        openPMD::Mesh mesh = iter.meshes[meshname];

        mesh.setGeometry(openPMD::Mesh::Geometry::cartesian);
        mesh.setAxisLabels(reversed(std::vector<std::string>{"x", "y", "z"}));
        mesh.setGridSpacing(to_vector<CCTK_REAL>(
            reversed(fmap([](auto x, auto y) { return x / CCTK_REAL(y); },
                          rdomain.hi - rdomain.lo, idomain.shape()))));
        mesh.setGridGlobalOffset(to_vector<double>(reversed(rdomain.lo)));
        mesh.setGridUnitSI(Unit::length);
        // const std::map<openPMD::UnitDimension, double> unitDimension{
        //     {openPMD::UnitDimension::L, 1}};
        // mesh.setUnitDimension(unitDimension);
        mesh.setTimeOffset(CCTK_REAL(0)); // TODO: check interface.ccl

        // Cell centred grids are offset by 1/2
        const Arith::vect<double, 3> position =
            fmap([](auto c) { return 0.5 * c; }, is_cell_centred);

        // Define tensor components

        // TODO: Set component names according to the tensor type
        std::vector<openPMD::MeshRecordComponent> record_components;
        record_components.reserve(numvars);
#if 0
        switch (numvars) {
        case 1:
          for (int vi = 0; vi < numvars; ++vi) {
           record_omponents.push_back(mesh[openPMD::MeshRecordComponent::SCALAR]);
          }
          break;
        case 3:
          for (int vi = 0; vi < numvars; ++vi) {
            const std::string cnames[] = {"x", "y", "z"};
            record_components.push_back(mesh[cnames[vi]]);
          }
          break;
        case 6:
          for (int vi = 0; vi < numvars; ++vi) {
            const std::string cnames[] = {"xx", "xy", "xz", "yy", "yz", "zz"};
            record_components.push_back(mesh[cnames[vi]]);
          }
          break;
        default:
          CCTK_VWARN(CCTK_WARN_ALERT, "unsupported tensor type gi=%d group=%s",
                     gi, CCTK_FullGroupName(gi));
          for (int vi = 0; vi < numvars; ++vi) {
          const std::string varname = CCTK_VarName(firstvarindex + vi);
            CCTK_VINFO("Creating component %d %s", vi, varname.c_str());
            record_components.push_back(mesh[varname]);
          }
          break;
        }
#endif
        for (int vi = 0; vi < numvars; ++vi) {
          const std::string componentname = make_componentname(gi, vi);
          record_components.push_back(mesh[componentname]);
          auto &record_component = record_components.back();
          record_component.setPosition(to_vector<double>(reversed(position)));
        }
        assert(int(record_components.size()) == numvars);

        // Write data

        if (io_verbose)
          CCTK_VINFO("Writing %d variables with %d components...", numvars,
                     num_local_components);

        for (int vi = 0; vi < numvars; ++vi)
          record_components.at(vi).resetDataset(dataset);

        // Loop over components (AMReX boxes)
        for (int local_component = 0; local_component < num_local_components;
             ++local_component) {
          const int component = mfab.IndexArray().at(local_component);

          const amrex::Box &fabbox =
              mfab.fabbox(component); // exterior (with ghosts)
          const box_t<int, 3> extbox{
            lo : {fabbox.smallEnd(0), fabbox.smallEnd(1), fabbox.smallEnd(2)},
            hi : {fabbox.bigEnd(0) + 1, fabbox.bigEnd(1) + 1,
                  fabbox.bigEnd(2) + 1}
          };
          const amrex::Box &validbox =
              mfab.box(component); // interior (without ghosts)
          const box_t<int, 3> intbox{
            lo : {validbox.smallEnd(0), validbox.smallEnd(1),
                  validbox.smallEnd(2)},
            hi : {validbox.bigEnd(0) + 1, validbox.bigEnd(1) + 1,
                  validbox.bigEnd(2) + 1}
          };
          const box_t<int, 3> &box = output_ghosts ? extbox : intbox;

          const openPMD::Offset start =
              to_vector(reversed(box.lo - idomain.lo));
          const openPMD::Extent count = to_vector(reversed(box.shape()));
          const int np = box.size();
          assert(int(count.at(0) * count.at(1) * count.at(2)) == np);
          for (int d = 0; d < 3; ++d)
            assert(start.at(d) >= 0);
          for (int d = 0; d < 3; ++d)
            assert(start.at(d) + count.at(d) <= extent.at(d));

          const amrex::FArrayBox &fab = mfab[component];
          for (int vi = 0; vi < numvars; ++vi) {
            if (output_ghosts) {
              const CCTK_REAL *const ptr = fab.dataPtr() + vi * np;
              record_components.at(vi).storeChunk(openPMD::shareRaw(ptr), start,
                                                  count);
            } else {
              std::shared_ptr<CCTK_REAL> ptr(
                  new CCTK_REAL[np], std::default_delete<CCTK_REAL[]>());
              const Arith::vect<int, 3> amrex_shape = extbox.shape();
              const Arith::vect<int, 3> amrex_offset = box.lo - extbox.lo;
              constexpr int amrex_di = 1;
              const int amrex_dj = amrex_di * amrex_shape[0];
              const int amrex_dk = amrex_dj * amrex_shape[1];
              const int amrex_np = amrex_dk * amrex_shape[2];
              const CCTK_REAL *restrict const amrex_ptr =
                  fab.dataPtr() + vi * amrex_np + amrex_di * amrex_offset[0] +
                  amrex_dj * amrex_offset[1] + amrex_dk * amrex_offset[2];
              const Arith::vect<int, 3> contig_shape = box.shape();
              constexpr int contig_di = 1;
              const int contig_dj = contig_di * contig_shape[0];
              const int contig_dk = contig_dj * contig_shape[1];
              const int contig_np = contig_dk * contig_shape[2];
              assert(contig_np == np);
              CCTK_REAL *restrict const contig_ptr = ptr.get();
              for (int k = 0; k < contig_shape[2]; ++k)
                for (int j = 0; j < contig_shape[1]; ++j)
                  for (int i = 0; i < contig_shape[0]; ++i)
                    contig_ptr[contig_di * i + contig_dj * j + contig_dk * k] =
                        amrex_ptr[amrex_di * i + amrex_dj * j + amrex_dk * k];
              record_components.at(vi).storeChunk(std::move(ptr), start, count);
            }
          } // for vi
        }   // for local_component
      }
    } // for gi

  } // for leveldata

  if (io_verbose)
    CCTK_VINFO("Closing iteration...");
  iter.close();

  if (CCTK_MyProc(nullptr) == 0) {
    std::ostringstream buf;
    buf << output_dir << "/" << output_file << ".openpmd.visit";
    const std::string visitname = buf.str();
    std::ofstream visit(visitname, ios::app);
    assert(visit.good());
    visit << *filename << "\n";
  }

  switch (iterationEncoding) {
  case openPMD::IterationEncoding::fileBased:
    write_iters.reset();
    series.reset();
    filename.reset();
    break;
  case openPMD::IterationEncoding::variableBased:
    // do nothing
    break;
  default:
    abort();
  }

  if (io_verbose)
    CCTK_VINFO("OutputOpenPMD done.");

  if (io_verbose)
    timer.print();
}

} // namespace CarpetX

#else

namespace CarpetX {
void ShutdownOpenPMD() {}
} // namespace CarpetX

#endif // #ifdef HAVE_CAPABILITY_OPENPMD
