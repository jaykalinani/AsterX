#ifndef EOS_UTILS_HXX
#define EOS_UTILS_HXX

#include <cmath>
#include <array>
#include <string>
#include <unordered_map>
#include <algorithm>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <loop_device.hxx>
#include <hdf5.h>
#include <mpi.h>

namespace EOSX {

using namespace std;
using namespace Loop;
using namespace Arith;

// Macro checking for errors coming from routines returning error codes
#define CHECK_ERROR(routine)                                                   \
  do {                                                                         \
    const auto err = routine;                                                  \
    if (err < 0) {                                                             \
      CCTK_VERROR("Routine '%s' returned error code %d", #routine, err);       \
    }                                                                          \
  } while (0)

/// Class representing a range
struct eos_range {
  CCTK_REAL min; ///< Minimum
  CCTK_REAL max; ///< Maximum
  eos_range() : min(0), max(0) {}
  eos_range(CCTK_REAL min_, CCTK_REAL max_) : min(min_), max(max_) {}
  bool contains(CCTK_REAL x) const { return (x >= min) && (x <= max); }
};

/// Class representing error conditions in EOS calls.
struct eos_status {
  bool failed; ///< Set to true if parameters are out of range/NAN/INF
  eos_status() : failed(false) {}
};

#ifdef H5_HAVE_PARALLEL
// Maps for warning messages (unchanged)
static const unordered_map<H5D_mpio_actual_io_mode_t,string> h5d_io_map = {
  {H5D_MPIO_NO_COLLECTIVE, "H5D_MPIO_NO_COLLECTIVE"},
  {H5D_MPIO_CHUNK_INDEPENDENT, "H5D_MPIO_CHUNK_INDEPENDENT"},
  {H5D_MPIO_CHUNK_COLLECTIVE, "H5D_MPIO_CHUNK_COLLECTIVE"},
  {H5D_MPIO_CHUNK_MIXED, "H5D_MPIO_CHUNK_MIXED"},
  {H5D_MPIO_CONTIGUOUS_COLLECTIVE, "H5D_MPIO_CONTIGUOUS_COLLECTIVE"}};
static const unordered_map<H5D_mpio_actual_chunk_opt_mode_t,string> h5d_chunk_map = {
  {H5D_MPIO_NO_CHUNK_OPTIMIZATION, "H5D_MPIO_NO_CHUNK_OPTIMIZATION"},
  {H5D_MPIO_LINK_CHUNK, "H5D_MPIO_LINK_CHUNK"},
  {H5D_MPIO_MULTI_CHUNK, "H5D_MPIO_MULTI_CHUNK"}};
static const unordered_map<uint32_t,string> h5p_cause_map = {
  {H5D_MPIO_COLLECTIVE, "H5D_MPIO_COLLECTIVE"},
  {H5D_MPIO_SET_INDEPENDENT, "H5D_MPIO_SET_INDEPENDENT"},
  {H5D_MPIO_DATATYPE_CONVERSION, "H5D_MPIO_DATATYPE_CONVERSION"},
  {H5D_MPIO_DATA_TRANSFORMS, "H5D_MPIO_DATA_TRANSFORMS"},
  {H5D_MPIO_MPI_OPT_TYPES_ENV_VAR_DISABLED, "H5D_MPIO_MPI_OPT_TYPES_ENV_VAR_DISABLED"},
  {H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES, "H5D_MPIO_NOT_SIMPLE_OR_SCALAR_DATASPACES"},
  {H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET, "H5D_MPIO_NOT_CONTIGUOUS_OR_CHUNKED_DATASET"},
  {H5D_MPIO_PARALLEL_FILTERED_WRITES_DISABLED, "H5D_MPIO_PARALLEL_FILTERED_WRITES_DISABLED"},
  {H5D_MPIO_ERROR_WHILE_CHECKING_COLLECTIVE_POSSIBLE, "H5D_MPIO_ERROR_WHILE_CHECKING_COLLECTIVE_POSSIBLE"},
  {H5D_MPIO_NO_COLLECTIVE_MAX_CAUSE, "H5D_MPIO_NO_COLLECTIVE_MAX_CAUSE"}};
#endif

// Read an integer dataset with *independent* MPI‐IO when in parallel:
inline void get_hdf5_int_dset(const hid_t &file_id,
                              const string &dset_name,
                              int npoints, int *var) {
  hid_t dset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
  assert(dset_id >= 0);

  hid_t dtype_id = H5Dget_type(dset_id);
  assert(dtype_id >= 0);
  assert(H5Tget_class(dtype_id) == H5T_INTEGER);
  H5Tclose(dtype_id);

  // create a transfer property list
  hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  assert(dxpl_id >= 0);
#ifdef H5_HAVE_PARALLEL
  // *** INDEPENDENT I/O ***
  CHECK_ERROR(H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT));
#endif

  CHECK_ERROR(H5Dread(dset_id, H5T_NATIVE_INT,
                      H5S_ALL, H5S_ALL,
                      dxpl_id, var));
  H5Pclose(dxpl_id);
  H5Dclose(dset_id);

#ifdef H5_HAVE_PARALLEL
  // (Optional) warn if still not truly independent
  H5D_mpio_actual_io_mode_t io_mode;
  CHECK_ERROR(H5Pget_mpio_actual_io_mode(dxpl_id, &io_mode));
  if (io_mode != H5D_MPIO_CHUNK_INDEPENDENT &&
      io_mode != H5D_MPIO_NO_COLLECTIVE) {
    CCTK_VWARN(1,
      "Dataset '%s' did not end up as independent I/O (mode=%s)",
      dset_name.c_str(), h5d_io_map.at(io_mode).c_str());
  }
#endif
}

// Read a real (double→CCTK_REAL) dataset with *independent* MPI‐IO when in parallel:
inline void get_hdf5_real_dset(const hid_t &file_id,
                               const string &dset_name,
                               int npoints, CCTK_REAL *var) {
  hid_t dset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
  assert(dset_id >= 0);

  // create a transfer property list
  hid_t dxpl_id = H5Pcreate(H5P_DATASET_XFER);
  assert(dxpl_id >= 0);
#ifdef H5_HAVE_PARALLEL
  // *** INDEPENDENT I/O ***
  CHECK_ERROR(H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_INDEPENDENT));
#endif

  CHECK_ERROR(H5Dread(dset_id, H5T_NATIVE_DOUBLE,
                      H5S_ALL, H5S_ALL,
                      dxpl_id, var));
  H5Pclose(dxpl_id);
  H5Dclose(dset_id);

#ifdef H5_HAVE_PARALLEL
  H5D_mpio_actual_io_mode_t io_mode;
  CHECK_ERROR(H5Pget_mpio_actual_io_mode(dxpl_id, &io_mode));
  if (io_mode != H5D_MPIO_CHUNK_INDEPENDENT &&
      io_mode != H5D_MPIO_NO_COLLECTIVE) {
    CCTK_VWARN(1,
      "Dataset '%s' did not end up as independent I/O (mode=%s)",
      dset_name.c_str(), h5d_io_map.at(io_mode).c_str());
  }
#endif
}

} // namespace EOSX

#endif // EOS_UTILS_HXX

