#ifndef EOS_3P_TABULATED3D_HXX
#define EOS_3P_TABULATED3D_HXX

#include <cctk.h>
#include <cmath>
#include <hdf5.h>
#include <mpi.h>
#include "eos_3p.hxx"
#include <string>
#include <AMReX.H>

using namespace std;

namespace EOSX {

class eos_3p_tabulated3d : public eos_3p {

//private:

  //linear_interp_uniform_ND_t<double, 1, 1> alltable;

public:
  CCTK_REAL gamma; // FIXME: get rid of this
  range rgeps;

  CCTK_INT ntemp, nrho, nye;
  // AMREX_GPU_MANAGED CCTK_REAL *logrho, *logtemp, *ye;

  // amrex::FArrayBox logtemp, logrho, ye;
  // amrex::Array4<CCTK_REAL> logpress, ...;

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  init(const range &rgeps_, const range &rgrho_, const range &rgye_) {
    set_range_rho(rgrho_);
    set_range_ye(rgye_);
    // TODO: first compute temp as a function of rho, ye, and eps, and then
    // initialize its range For now, as dummy, we pass range of eps as range of
    // temp
    set_range_temp(rgeps_);
  }

  // Routine reading an HDF5 simple dataset consisting of one integer element
  // FIXME: make this more general: any type (template), any number of elements
  // (but still 1D-shaped)
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_INT
  get_hdf5_simple_int(const hid_t &file_id, const string &dset_name) {
    const auto dset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    assert(dset_id >= 0);

    const auto dspace_id = H5Dget_space(dset_id);
    assert(dspace_id >= 0);

    const auto dspacetype_id = H5Sget_simple_extent_type(dspace_id);
    assert(dspacetype_id == H5S_SIMPLE);

    const auto ndims = H5Sget_simple_extent_ndims(dspace_id);
    assert(ndims == 1);

    hsize_t size;
    const auto ndims_again =
        H5Sget_simple_extent_dims(dspace_id, &size, nullptr);
    CHECK_ERROR(H5Sclose(dspace_id));
    assert(ndims_again == 1);
    assert(size == 1);

    const auto dtype_id = H5Dget_type(dset_id);
    assert(dtype_id >= 0);

    const auto dtypeclass = H5Tget_class(dtype_id);
    assert(dtypeclass == H5T_INTEGER);
    CHECK_ERROR(H5Tclose(dtype_id));

    auto dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    assert(dxpl_id >= 0);

#ifdef H5_HAVE_PARALLEL
    CHECK_ERROR(H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE));
#else
    dxpl_id = H5P_DEFAULT;
#endif

    CCTK_INT buffer;
    CHECK_ERROR(
        H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, dxpl_id, &buffer));
    CHECK_ERROR(H5Dclose(dset_id));

#ifdef H5_HAVE_PARALLEL
    H5D_mpio_actual_io_mode_t actual_io_mode;
    CHECK_ERROR(H5Pget_mpio_actual_io_mode(dxpl_id, &actual_io_mode));

    H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
    CHECK_ERROR(
        H5Pget_mpio_actual_chunk_opt_mode(dxpl_id, &actual_chunk_opt_mode));

    uint32_t no_collective_cause_local, no_collective_cause_global;
    CHECK_ERROR(H5Pget_mpio_no_collective_cause(
        dxpl_id, &no_collective_cause_local, &no_collective_cause_global));

    if (actual_io_mode == H5D_MPIO_NO_COLLECTIVE or
        // actual_chunk_opt_mode      == H5D_MPIO_NO_CHUNK_OPTIMIZATION or  //
        // In general, input files are not chunked
        no_collective_cause_local != H5D_MPIO_COLLECTIVE or
        no_collective_cause_global != H5D_MPIO_COLLECTIVE) {
      CCTK_VWARN(
          1,
          "Actual I/O mode, chunk optimization and local and global "
          "non-collective I/O causes when reading data from dataset '%s': "
          "'%s', '%s', '%s', '%s'",
          H5D_mpio_actual_io_mode_map.at(actual_io_mode).c_str(),
          H5D_mpio_actual_chunk_opt_mode_map.at(actual_chunk_opt_mode).c_str(),
          H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_local)
              .c_str(),
          H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_global)
              .c_str(),
          dset_name.c_str());
    }
#endif // H5_HAVE_PARALLEL

    CHECK_ERROR(H5Pclose(dxpl_id));
    return buffer;
  }

  // Routine reading a 1D HDF5 simple dataset storing real numbers
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  get_hdf5_simple_1Darray(const hid_t &file_id, const string &dset_name,
                          amrex::FArrayBox &var) {
    const auto dset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    assert(dset_id >= 0);

    const auto dspace_id = H5Dget_space(dset_id);
    assert(dspace_id >= 0);

    const auto dspacetype_id = H5Sget_simple_extent_type(dspace_id);
    assert(dspacetype_id == H5S_SIMPLE);

    const auto ndims = H5Sget_simple_extent_ndims(dspace_id);
    assert(ndims == 1);

    hsize_t size;
    const auto ndims_again =
        H5Sget_simple_extent_dims(dspace_id, &size, nullptr);
    CHECK_ERROR(H5Sclose(dspace_id));
    assert(ndims_again == 1);
    assert(size > 0);

    const auto dtype_id = H5Dget_type(dset_id);
    assert(dtype_id >= 0);

    const auto dtypeclass = H5Tget_class(dtype_id);
    assert(dtypeclass == H5T_FLOAT);
    CHECK_ERROR(H5Tclose(dtype_id));

    auto dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    assert(dxpl_id >= 0);

#ifdef H5_HAVE_PARALLEL
    CHECK_ERROR(H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE));
#else
    dxpl_id = H5P_DEFAULT;
#endif

    CCTK_REAL buffer[size];
    // CCTK_REAL *buffer = new CCTK_REAL[size];
    CHECK_ERROR(
        H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl_id, buffer));
    CHECK_ERROR(H5Dclose(dset_id));

#ifdef H5_HAVE_PARALLEL
    H5D_mpio_actual_io_mode_t actual_io_mode;
    CHECK_ERROR(H5Pget_mpio_actual_io_mode(dxpl_id, &actual_io_mode));

    H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
    CHECK_ERROR(
        H5Pget_mpio_actual_chunk_opt_mode(dxpl_id, &actual_chunk_opt_mode));

    uint32_t no_collective_cause_local, no_collective_cause_global;
    CHECK_ERROR(H5Pget_mpio_no_collective_cause(
        dxpl_id, &no_collective_cause_local, &no_collective_cause_global));

    if (actual_io_mode == H5D_MPIO_NO_COLLECTIVE or
        // actual_chunk_opt_mode      == H5D_MPIO_NO_CHUNK_OPTIMIZATION or  //
        // In general, input files are not chunked
        no_collective_cause_local != H5D_MPIO_COLLECTIVE or
        no_collective_cause_global != H5D_MPIO_COLLECTIVE) {
      CCTK_VWARN(
          1,
          "Actual I/O mode, chunk optimization and local and global "
          "non-collective I/O causes when reading data from dataset '%s': "
          "'%s', '%s', '%s', '%s'",
          H5D_mpio_actual_io_mode_map.at(actual_io_mode).c_str(),
          H5D_mpio_actual_chunk_opt_mode_map.at(actual_chunk_opt_mode).c_str(),
          H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_local)
              .c_str(),
          H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_global)
              .c_str(),
          dset_name.c_str());
    }
#endif // H5_HAVE_PARALLEL

    CHECK_ERROR(H5Pclose(dxpl_id));

    // amrex::Box box1d(amrex::IntVect{0, 0, 0}, amrex::IntVect{size - 1, 0,
    // 0}); var = amrex::FArrayBox(box1d, 1, amrex::The_Managed_Arena());
    // TODO: fill var with buffer

    return;
  }

  // Routine reading the EOS table and filling the corresponding object
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  read_eos_table(const string &filename) {
    CCTK_VINFO("Reading EOS table '%s'", filename.c_str());

    auto fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    assert(fapl_id >= 0);
    hid_t file_id = 0;
    int rank_id;
    CHECK_ERROR(MPI_Comm_rank(MPI_COMM_WORLD, &rank_id));

#ifdef H5_HAVE_PARALLEL
    CHECK_ERROR(H5Pset_fapl_mpio(fapl_id, MPI_COMM_WORLD, MPI_INFO_NULL));
    CHECK_ERROR(H5Pset_all_coll_metadata_ops(fapl_id, true));
    file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id);
#else
    fapl_id = H5P_DEFAULT;
    if (rank_id == 0) {
      file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, fapl_id);
    }
// TODO: bcast table info to all ranks
#endif

    assert(file_id >= 0);

    // TODO: finish reading the EOS table and filling the EOS object
    ntemp = get_hdf5_simple_int(file_id, "pointstemp");
    nrho = get_hdf5_simple_int(file_id, "pointsrho");
    nye = get_hdf5_simple_int(file_id, "pointsye");

    CCTK_VINFO("EOS table dimensions: ntemp = %d, nrho = %d, nye = %d", ntemp,
               nrho, nye);

    CHECK_ERROR(H5Pclose(fapl_id));

#ifdef H5_HAVE_PARALLEL
    CHECK_ERROR(H5Fclose(file_id));
#else
    if (rank_id == 0) {
      CHECK_ERROR(H5Fclose(file_id));
    }
#endif

    return;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                              const CCTK_REAL ye) const {
    CCTK_REAL press = 0.0; // tab3d_press(rho, eps, ye, &ierr);

    //
    //press = alltable.interpolate<PRESS_>(rho,T,ye)[0]
    //

    return press;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_valid_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                             const CCTK_REAL ye) const {
    CCTK_REAL csnd2 = 0.0; // tab3d_csnd2(rho, eps, ye, &ierr);
    assert(csnd2 >= 0);    // Soundspeed^2 should never ever be negative
    return sqrt(csnd2);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  temp_from_valid_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps,
                             const CCTK_REAL ye) const {
    CCTK_REAL temp = 0.0; // tab3d_temp(rho, eps, ye, &ierr);
    return temp;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void
  press_derivs_from_valid_rho_eps_ye(CCTK_REAL &press, CCTK_REAL &dpdrho,
                                     CCTK_REAL &dpdeps, const CCTK_REAL rho,
                                     const CCTK_REAL eps,
                                     const CCTK_REAL ye) const {
    printf("press_derivs_from_valid_rho_eps_ye is not supported anymore!");
    exit(EXIT_FAILURE);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  press_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                               const CCTK_REAL ye) const {
    return 0.0; // tab3d_press_from_temp(rho, temp, ye);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  csnd_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                              const CCTK_REAL ye) const {
    CCTK_REAL csnd2 = 0.0; // tab3d_csnd2_from_temp(rho, temp, ye);
    assert(csnd2 >= 0);    // Soundspeed^2 should never ever be negative
    return sqrt(csnd2);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  entropy_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp,
                                 const CCTK_REAL ye) const {
    return 0.0; // tab3d_entropy_from_temp(rho, temp, ye);
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline range
  range_eps_from_valid_rho_ye(const CCTK_REAL rho, const CCTK_REAL ye) const {
    return rgeps;
  }

  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
  eps_from_valid_rho_press_ye(const CCTK_REAL rho, const CCTK_REAL press,
                              const CCTK_REAL ye) const {
    return 0.0; // press / (rho * gm1);
  }
};
} // namespace EOSX

#endif
