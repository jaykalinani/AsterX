#ifndef EOS_3P_TABULATED3D_HXX
#define EOS_3P_TABULATED3D_HXX

#define NTABLES 19
#define LENGTHGF 6.77269222552442e-06
#define TIMEGF 2.03040204956746e05
#define RHOGF 1.61887093132742e-18
#define PRESSGF 1.80123683248503e-39
#define EPSGF 1.11265005605362e-21
#define INVRHOGF 6.17714470405638e17
#define INVEPSGF 8.98755178736818e20
#define INVPRESSGF 5.55174079257738e38

#include <cctk.h>
#include <cmath>
#include <hdf5.h>
#include <mpi.h>
#include "eos_3p.hxx"
#include <string>
#include <AMReX.H>

using namespace std;

namespace EOSX {

using namespace amrex;

class eos_3p_tabulated3d : public eos_3p {

public:
  CCTK_REAL gamma;  // FIXME: get rid of this
  range rgeps;
  
  CCTK_INT ntemp, nrho, nye;

  CCTK_REAL *logrho, *logtemp, *yes; // FIXME: AMREX_GPU_MANAGED?
  CCTK_REAL *alltables;
  CCTK_REAL *epstable;
  CCTK_REAL energy_shift;
  
  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void init(
    const range &rgeps_, const range &rgrho_, const range &rgye_)
  {
    set_range_rho(rgrho_);
    set_range_ye(rgye_);
    // TODO: first compute temp as a function of rho, ye, and eps, and then initialize its range
    // For now, as dummy, we pass range of eps as range of temp
    set_range_temp(rgeps_);
  }


  // Routine reading an HDF5 integer dataset
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline
  void get_hdf5_int_dset(const hid_t  &file_id,
             const string &dset_name,
			       const int npoints,
			       int* var) {

    const auto dset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    assert(dset_id >= 0);

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

    CHECK_ERROR(H5Dread(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, dxpl_id, var));
    CHECK_ERROR(H5Dclose(dset_id));

    #ifdef H5_HAVE_PARALLEL
    H5D_mpio_actual_io_mode_t actual_io_mode;
    CHECK_ERROR(H5Pget_mpio_actual_io_mode(dxpl_id, &actual_io_mode));

    H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
    CHECK_ERROR(H5Pget_mpio_actual_chunk_opt_mode(dxpl_id, &actual_chunk_opt_mode));

    uint32_t no_collective_cause_local, no_collective_cause_global;
    CHECK_ERROR(H5Pget_mpio_no_collective_cause(dxpl_id,
                &no_collective_cause_local, &no_collective_cause_global));

    if (actual_io_mode             == H5D_MPIO_NO_COLLECTIVE         or
        //actual_chunk_opt_mode      == H5D_MPIO_NO_CHUNK_OPTIMIZATION or  // In general, input files are not chunked
        no_collective_cause_local  != H5D_MPIO_COLLECTIVE            or
        no_collective_cause_global != H5D_MPIO_COLLECTIVE)
    {
        CCTK_VWARN(1, "Actual I/O mode, chunk optimization and local and global non-collective I/O causes when reading data from dataset '%s': '%s', '%s', '%s', '%s'",
                   H5D_mpio_actual_io_mode_map.at(actual_io_mode).c_str(),
                   H5D_mpio_actual_chunk_opt_mode_map.at(actual_chunk_opt_mode).c_str(),
                   H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_local).c_str(),
                   H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_global).c_str(),
                   dset_name.c_str());
    }
    #endif  // H5_HAVE_PARALLEL

    CHECK_ERROR(H5Pclose(dxpl_id));
  }

  // Routine reading an HDF5 real number dataset
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline
  void get_hdf5_real_dset(const hid_t  &file_id,
             const string &dset_name,
			       const int npoints,
			       double* var) {

    const auto dset_id = H5Dopen(file_id, dset_name.c_str(), H5P_DEFAULT);
    assert(dset_id >= 0);

    auto dxpl_id = H5Pcreate(H5P_DATASET_XFER);
    assert(dxpl_id >= 0);

    #ifdef H5_HAVE_PARALLEL
    CHECK_ERROR(H5Pset_dxpl_mpio(dxpl_id, H5FD_MPIO_COLLECTIVE));
    #else
    dxpl_id = H5P_DEFAULT;
    #endif

    CHECK_ERROR(H5Dread(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, dxpl_id, var));
    CHECK_ERROR(H5Dclose(dset_id));

    #ifdef H5_HAVE_PARALLEL
    H5D_mpio_actual_io_mode_t actual_io_mode;
    CHECK_ERROR(H5Pget_mpio_actual_io_mode(dxpl_id, &actual_io_mode));

    H5D_mpio_actual_chunk_opt_mode_t actual_chunk_opt_mode;
    CHECK_ERROR(H5Pget_mpio_actual_chunk_opt_mode(dxpl_id, &actual_chunk_opt_mode));

    uint32_t no_collective_cause_local, no_collective_cause_global;
    CHECK_ERROR(H5Pget_mpio_no_collective_cause(dxpl_id,
                &no_collective_cause_local, &no_collective_cause_global));

    if (actual_io_mode             == H5D_MPIO_NO_COLLECTIVE         or
        //actual_chunk_opt_mode      == H5D_MPIO_NO_CHUNK_OPTIMIZATION or  // In general, input files are not chunked
        no_collective_cause_local  != H5D_MPIO_COLLECTIVE            or
        no_collective_cause_global != H5D_MPIO_COLLECTIVE)
    {
        CCTK_VWARN(1, "Actual I/O mode, chunk optimization and local and global non-collective I/O causes when reading data from dataset '%s': '%s', '%s', '%s', '%s'",
                   H5D_mpio_actual_io_mode_map.at(actual_io_mode).c_str(),
                   H5D_mpio_actual_chunk_opt_mode_map.at(actual_chunk_opt_mode).c_str(),
                   H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_local).c_str(),
                   H5Pget_mpio_no_collective_cause_map.at(no_collective_cause_global).c_str(),
                   dset_name.c_str());
    }
    #endif  // H5_HAVE_PARALLEL

    CHECK_ERROR(H5Pclose(dxpl_id));
  }


  // Routine reading the EOS table and filling the corresponding object
  CCTK_HOST CCTK_ATTRIBUTE_ALWAYS_INLINE inline void read_eos_table(const string &filename) {
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

		// Get number of points
    get_hdf5_int_dset(file_id, "pointstemp", 1, &ntemp);
    get_hdf5_int_dset(file_id, "pointsrho", 1, &nrho);
    get_hdf5_int_dset(file_id, "pointsye", 1, &nye);
 
    const int npoints = ntemp * nrho * nye;

    CCTK_VINFO("EOS table dimensions: ntemp = %d, nrho = %d, nye = %d", ntemp, nrho, nye);

    // Allocate memory for tables
    double* alltables_temp;
    if (!(alltables_temp = (double*)The_Managed_Arena()->alloc(npoints * NTABLES * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
    	"Cannot allocate memory for EOS table");
    }
    if (!(logrho = (double*)The_Managed_Arena()->alloc(nrho * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
    	"Cannot allocate memory for EOS table");
    }
    if (!(logtemp = (double*)The_Managed_Arena()->alloc(ntemp * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
    	"Cannot allocate memory for EOS table");
    }
    if (!(yes = (double*)The_Managed_Arena()->alloc(nye * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
    	"Cannot allocate memory for EOS table");
    }

    // Prepare HDF5 to read hyperslabs into alltables_temp
    hsize_t table_dims[2] = {NTABLES, (hsize_t)npoints};
    hsize_t var3[2]       = { 1, (hsize_t)npoints};
    hid_t mem3 =  H5Screate_simple(2, table_dims, NULL);

    // hydro (and munu)
    get_hdf5_real_dset(file_id, "logpress",  npoints, &alltables_temp[0  * npoints]);
    get_hdf5_real_dset(file_id, "logenergy", npoints, &alltables_temp[1  * npoints]);
    get_hdf5_real_dset(file_id, "entropy",   npoints, &alltables_temp[2  * npoints]);
    get_hdf5_real_dset(file_id, "munu",      npoints, &alltables_temp[3  * npoints]);
    get_hdf5_real_dset(file_id, "cs2",       npoints, &alltables_temp[4  * npoints]);
    get_hdf5_real_dset(file_id, "dedt",      npoints, &alltables_temp[5  * npoints]);
    get_hdf5_real_dset(file_id, "dpdrhoe",   npoints, &alltables_temp[6  * npoints]);
    get_hdf5_real_dset(file_id, "dpderho",   npoints, &alltables_temp[7  * npoints]);

    // chemical potentials
    get_hdf5_real_dset(file_id, "muhat",     npoints, &alltables_temp[8  * npoints]);
    get_hdf5_real_dset(file_id, "mu_e",      npoints, &alltables_temp[9  * npoints]);
    get_hdf5_real_dset(file_id, "mu_p",			 npoints, &alltables_temp[10 * npoints]);
    get_hdf5_real_dset(file_id, "mu_n",			 npoints, &alltables_temp[11 * npoints]);
    
    // compositions
    get_hdf5_real_dset(file_id, "Xa", 			 npoints, &alltables_temp[12 * npoints]);
    get_hdf5_real_dset(file_id, "Xh", 			 npoints, &alltables_temp[13 * npoints]);
    get_hdf5_real_dset(file_id, "Xn", 			 npoints, &alltables_temp[14 * npoints]);
    get_hdf5_real_dset(file_id, "Xp", 			 npoints, &alltables_temp[15 * npoints]);

    // average nucleus
    get_hdf5_real_dset(file_id, "Abar",      npoints, &alltables_temp[16 * npoints]);
    get_hdf5_real_dset(file_id, "Zbar",      npoints, &alltables_temp[17 * npoints]);

    // Gamma
    get_hdf5_real_dset(file_id, "gamma",		 npoints, &alltables_temp[18 * npoints]);

    // Read additional tables and variables
    get_hdf5_real_dset(file_id, "logrho",		  nrho, logrho);
    get_hdf5_real_dset(file_id, "logtemp",	  ntemp, logtemp);
    get_hdf5_real_dset(file_id, "ye",				  nye, yes);
    get_hdf5_real_dset(file_id, "energy_shift", 1, &energy_shift);

    CHECK_ERROR(H5Pclose(fapl_id));
    CHECK_ERROR(H5Sclose(mem3));

    #ifdef H5_HAVE_PARALLEL
    CHECK_ERROR(H5Fclose(file_id));
    #else
    if (rank_id == 0) {
      CHECK_ERROR(H5Fclose(file_id));
    }
    #endif

    // Fill actual table
    if (!(alltables = (double*)The_Managed_Arena()->alloc(npoints * NTABLES * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
    	"Cannot allocate memory for EOS table");
    }
    for(int iv = 0;iv<NTABLES;iv++)
    	for(int k = 0; k<nye;k++)
      	    for(int j = 0; j<ntemp; j++)
  		for(int i = 0; i<nrho; i++) {
	    int indold = i + nrho*(j + ntemp*(k + nye*iv));
	    int indnew = iv + NTABLES*(i + nrho*(j + ntemp*k));

	    // Maybe swap temp axis?
	    // int indnew = iv + NTABLES*(j + ntemp*(i + nrho*k));
	    alltables[indnew] = alltables_temp[indold];
    }

    // free memory of temporary array
    The_Managed_Arena()->free(alltables_temp);

  	// allocate epstable; a linear-scale eps table
    // that allows us to extrapolate to negative eps
    if (!(epstable = (double*)The_Managed_Arena()->alloc(npoints * sizeof(double)))) {
    CCTK_VError(__LINE__, __FILE__, CCTK_THORNSTRING,
    	"Cannot allocate memory for EOS table");
    }
  
    // convert units, convert logs to natural log
    // The latter is great, because exp() is way faster than pow()
    // pressure
    energy_shift = energy_shift * EPSGF;
    for(int i=0;i<nrho;i++) {
      // rewrite:
      //logrho[i] = log(pow(10.0,logrho[i]) * RHOGF);
      // by using log(a^b*c) = b*log(a)+log(c)
      logrho[i] = logrho[i] * log(10.) + log(RHOGF);
    }

    for(int i=0;i<ntemp;i++) {
      //logtemp[i] = log(pow(10.0,logtemp[i]));
      logtemp[i] = logtemp[i]*log(10.0);
    }

      // convert units
    for(int i=0;i<npoints;i++) {

      { // pressure
        int idx = 0 + NTABLES*i;
        alltables[idx] = alltables[idx] * log(10.0) + log(PRESSGF);
      }

      { // eps
        int idx = 1 + NTABLES*i;
        alltables[idx] = alltables[idx] * log(10.0) + log(EPSGF);
        epstable[i] = exp(alltables[idx]);
      }

      { // cs2
        int idx = 4 + NTABLES*i;
        alltables[idx] *= LENGTHGF*LENGTHGF/TIMEGF/TIMEGF;
      }

      { // dedT
        int idx = 5 + NTABLES*i;
        alltables[idx] *= EPSGF;
      }

      { // dpdrhoe
        int idx = 6 + NTABLES*i;
        alltables[idx] *= PRESSGF/RHOGF;
      }

      { // dpderho
        int idx = 7 + NTABLES*i;
        alltables[idx] *= PRESSGF/EPSGF;
      }

    }

    // set up steps, mins, maxes here?
    return;
  }


  CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL press_from_valid_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps, const CCTK_REAL ye) const
{
  CCTK_REAL press = 0.0; //tab3d_press(rho, eps, ye, &ierr);
  return press;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL csnd_from_valid_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps, const CCTK_REAL ye) const
{
  CCTK_REAL csnd2 = 0.0; //tab3d_csnd2(rho, eps, ye, &ierr);
  assert(csnd2 >= 0); //Soundspeed^2 should never ever be negative
  return sqrt(csnd2);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL temp_from_valid_rho_eps_ye(const CCTK_REAL rho, const CCTK_REAL eps, const CCTK_REAL ye) const
{
  CCTK_REAL temp = 0.0; //tab3d_temp(rho, eps, ye, &ierr);
  return temp;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline void press_derivs_from_valid_rho_eps_ye(CCTK_REAL& press, CCTK_REAL& dpdrho, CCTK_REAL& dpdeps,
           const CCTK_REAL rho, const CCTK_REAL eps, const CCTK_REAL ye) const
{
  printf("press_derivs_from_valid_rho_eps_ye is not supported anymore!");
  exit(EXIT_FAILURE);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL press_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp, const CCTK_REAL ye) const
{
    return 0.0; //tab3d_press_from_temp(rho, temp, ye);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL csnd_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp, const CCTK_REAL ye) const
{
  CCTK_REAL csnd2 = 0.0; //tab3d_csnd2_from_temp(rho, temp, ye);
  assert(csnd2 >= 0); //Soundspeed^2 should never ever be negative
  return sqrt(csnd2);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL entropy_from_valid_rho_temp_ye(const CCTK_REAL rho, const CCTK_REAL temp, const CCTK_REAL ye) const
{
    return 0.0; //tab3d_entropy_from_temp(rho, temp, ye);
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline range
 range_eps_from_valid_rho_ye(const CCTK_REAL rho,
                                          const CCTK_REAL ye) const {
  return rgeps;
}

CCTK_HOST CCTK_DEVICE CCTK_ATTRIBUTE_ALWAYS_INLINE inline CCTK_REAL
eps_from_valid_rho_press_ye(const CCTK_REAL rho,
                                          const CCTK_REAL press,
                                          const CCTK_REAL ye) const {
  return 0.0; //press / (rho * gm1);
}

};
} // namespace EOSX

#endif
