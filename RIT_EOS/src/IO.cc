#include <cmath>      // For log() and exp()
#include <H5Cpp.h>
#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include "RIT_EOS.hh"

#ifdef HAVE_CAPABILITY_MPI
#include <mpi.h>
#endif

// ***** DEBUG ONLY *****
//#include <iostream>
//#include <vector>



// Definition of namespace "RIT_EOS_tablevars" (declared in RIT_EOS.hh)
namespace RIT_EOS_tablevars {
    CCTK_INT  nrho,    ntemp,    nye;
    CCTK_INT  nrho_m1, ntemp_m1, nye_m1;
    CCTK_REAL energy_shift;
    CCTK_REAL *restrict logrho,   *restrict logtemp, *restrict ye;
    CCTK_REAL *restrict logpress, *restrict logenergy;
    CCTK_REAL *restrict dpdrhoe,  *restrict dpderho;
}

using H5::H5File;
using H5::DataSet;
using H5::PredType;
using namespace RIT_EOS_tablevars;
using RIT_EOS_helpers::ln10;
using RIT_EOS_helpers::eps_CGStoGeom;
using RIT_EOS_helpers::log_rho_CGStoGeom;
using RIT_EOS_helpers::log_press_CGStoGeom;
using RIT_EOS_helpers::log_eps_CGStoGeom;
using RIT_EOS_helpers::press_over_rho_CGStoGeom;
using RIT_EOS_helpers::press_over_eps_CGStoGeom;

using namespace RIT_EOS_keys;
key RIT_EOS_keys::press_key;
key RIT_EOS_keys::eps_key;
key RIT_EOS_keys::dpdrhoe_key;
key RIT_EOS_keys::dpderho_key;

// ***** DEBUG ONLY *****
//using std::vector;





// Helper routine to get nrho, ntemp, nye from the EOS table
void read_nrho_ntemp_nye_from_EOStable(H5File &file) {
    // Declare the DataSet variable
    DataSet dataset;

    // Get nrho, ntemp, nye from the EOS table
    dataset = file.openDataSet("pointsrho");
    dataset.read(&nrho, PredType::NATIVE_INT, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'pointsrho'    successfully read from the EOS table");

    dataset = file.openDataSet("pointstemp");
    dataset.read(&ntemp, PredType::NATIVE_INT, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'pointstemp'   successfully read from the EOS table");

    dataset = file.openDataSet("pointsye");
    dataset.read(&nye, PredType::NATIVE_INT, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'pointsye'     successfully read from the EOS table");


    // Sanity checks
    if (nrho  < 1) CCTK_VERROR("nrho = %d < 1",  nrho);
    if (ntemp < 1) CCTK_VERROR("ntemp = %d < 1", ntemp);
    if (nye   < 1) CCTK_VERROR("nye = %d < 1",   nye);

    // Set helper variables
    nrho_m1  = nrho  - 1;
    ntemp_m1 = ntemp - 1;
    nye_m1   = nye   - 1;

    // Print info about the dimensions of the EOS table
    CCTK_VINFO("EOS table dimensions: nrho = %d, ntemp = %d, nye = %d",
               nrho, ntemp, nye);

    return;
}





// Helper routine to get the array variables from the EOS table
void read_arrays_from_EOStable(H5File &file) {
    // Declare the DataSet variable
    DataSet dataset;

    // Get array quantities from the EOS table
    dataset = file.openDataSet("logrho");
    dataset.read(logrho, PredType::NATIVE_DOUBLE, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'logrho'       successfully read from the EOS table");

    dataset = file.openDataSet("logtemp");
    dataset.read(logtemp, PredType::NATIVE_DOUBLE, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'logtemp'      successfully read from the EOS table");

    dataset = file.openDataSet("ye");
    dataset.read(ye, PredType::NATIVE_DOUBLE, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'ye'           successfully read from the EOS table");

    dataset = file.openDataSet("logpress");
    dataset.read(logpress, PredType::NATIVE_DOUBLE, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'logpress'     successfully read from the EOS table");

    dataset = file.openDataSet("logenergy");
    dataset.read(logenergy, PredType::NATIVE_DOUBLE, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'logenergy'    successfully read from the EOS table");

    dataset = file.openDataSet("dpdrhoe");
    dataset.read(dpdrhoe, PredType::NATIVE_DOUBLE, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'dpdrhoe'      successfully read from the EOS table");

    dataset = file.openDataSet("dpderho");
    dataset.read(dpderho, PredType::NATIVE_DOUBLE, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'dpderho'      successfully read from the EOS table");

    dataset = file.openDataSet("energy_shift");
    dataset.read(&energy_shift, PredType::NATIVE_DOUBLE, H5S_ALL, H5S_ALL);
    CCTK_INFO("Dataset 'energy_shift' successfully read from the EOS table");

    return;
}





// Routine to read the EOS table
extern "C" void ReadEOStable(CCTK_ARGUMENTS) {
    DECLARE_CCTK_PARAMETERS;

    // Declaration of variables
    CCTK_INT  i, j, k, kji;
    CCTK_INT  procID, Bcast_successful;
    CCTK_INT  nrho_ntemp_nye, nrho_ntemp_nye_buffer[3];
    H5File    file;

    // Print some info
    CCTK_INFO("The following conversion factors between EOS table units (mostly CGS) and geometrized (\"Cactus\") units will be used:");
    CCTK_VINFO("    Mass density:             %.6e (approximately)", exp(log_rho_CGStoGeom));
    CCTK_VINFO("    Pressure:                 %.6e (approximately)", exp(log_press_CGStoGeom));
    CCTK_VINFO("    Specific internal energy: %.6e (approximately)", eps_CGStoGeom);


    /* FIXME: if I don't specify any '#pragma omp' directive, will the following
              code be executed by multiple threads? I don't want this, because
              I/O operations should only be executed by one single thread!      */


    #ifdef HAVE_CAPABILITY_MPI
        /* If MPI is enabled, then only do I/O operations on the master MPI
           process and later broadcast the data to all other MPI processes with
           MPI_Bcast()                                                          */
        MPI_Comm_rank(MPI_COMM_WORLD, &procID);  // Get MPI process ID

        if (procID == 0) {
            // Open HDF5 EOS table file
            file.openFile(EOStable_path, H5F_ACC_RDONLY, H5P_DEFAULT);

            // Get nrho, ntemp, nye from the EOS table
            read_nrho_ntemp_nye_from_EOStable(file);

            /* Fill a buffer array with nrho, ntemp, nye. This buffer is used
               below to broadcast nrho, ntemp, nye from the master MPI process
               to all other MPI processes.                                      */
            nrho_ntemp_nye_buffer[0] = nrho;
            nrho_ntemp_nye_buffer[1] = ntemp;
            nrho_ntemp_nye_buffer[2] = nye;
        }

        // Synchronize all MPI processes // FIXME: is this really needed?
        MPI_Barrier(MPI_COMM_WORLD);

        /* Broadcast nrho, ntemp, nye from the master process (0) to all other
           MPI processes                                                        */
        Bcast_successful = MPI_Bcast(nrho_ntemp_nye_buffer, 3, MPI_INT, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast nrho, ntemp, nye from the master MPI process to all other MPI processes");

        // Synchronize all MPI processes // FIXME: is this really needed?
        MPI_Barrier(MPI_COMM_WORLD);

        /* Now all non-master MPI processes must retrieve nrho, ntemp, nye from
           buffer 'nrho_ntemp_nye_buffer'.                                      */
        nrho  = nrho_ntemp_nye_buffer[0];
        ntemp = nrho_ntemp_nye_buffer[1];
        nye   = nrho_ntemp_nye_buffer[2];

        // Synchronize all MPI processes // FIXME: is this really needed?
        MPI_Barrier(MPI_COMM_WORLD);

    #else
        // Open HDF5 EOS table file
        file.openFile(EOStable_path, H5F_ACC_RDONLY, H5P_DEFAULT);

        // Get nrho, ntemp, nye from the EOS table
        read_nrho_ntemp_nye_from_EOStable(file);
    #endif



    /* Allocate memory for the quantities to be read from the table
       N.B.: memory must be allocated on ALL MPI processes, not just on the
             master process!                                                    */
    nrho_ntemp_nye = nye*ntemp*nrho;

    logrho    = new double[nrho];
    logtemp   = new double[ntemp];
    ye        = new double[nye];
    logpress  = new double[nrho_ntemp_nye];
    logenergy = new double[nrho_ntemp_nye];
    dpdrhoe   = new double[nrho_ntemp_nye];
    dpderho   = new double[nrho_ntemp_nye];

    /* Check that the memory for the quantities to be read from the EOS
       table has been allocated correctly                                       */
    /* FIXME: checking if a pointer is NULL does NOT seem to catch failures in
              memory allocation. How can I catch those errors?                  */
    if (logrho    == nullptr)  CCTK_ERROR("Failed to allocate memory for pointer 'logrho'");
    if (logtemp   == nullptr)  CCTK_ERROR("Failed to allocate memory for pointer 'logtemp'");
    if (ye        == nullptr)  CCTK_ERROR("Failed to allocate memory for pointer 'ye'");
    if (logpress  == nullptr)  CCTK_ERROR("Failed to allocate memory for pointer 'logpress'");
    if (logenergy == nullptr)  CCTK_ERROR("Failed to allocate memory for pointer 'logenergy'");
    if (dpdrhoe   == nullptr)  CCTK_ERROR("Failed to allocate memory for pointer 'dpdrhoe'");
    if (dpderho   == nullptr)  CCTK_ERROR("Failed to allocate memory for pointer 'dpderho'");



    #ifdef HAVE_CAPABILITY_MPI
        // Synchronize all MPI processes
        /* FIXME: is this really needed? And in any case, there's a barrier at
                  end of the previous MPI if statement...                       */
        MPI_Barrier(MPI_COMM_WORLD);

        // Get array quantities from the EOS table
        if (procID == 0)
            read_arrays_from_EOStable(file);

        // Synchronize all MPI processes // FIXME: is this really needed?
        MPI_Barrier(MPI_COMM_WORLD);


        /* Broadcast all arrays from the master process (0) to all other MPI
           processes                                                            */
        Bcast_successful = MPI_Bcast(logrho, nrho, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast array 'logrho' from the master MPI process to all other MPI processes");

        Bcast_successful = MPI_Bcast(logtemp, ntemp, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast array 'logtemp' from the master MPI process to all other MPI processes");

        Bcast_successful = MPI_Bcast(ye, nye, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast array 'ye' from the master MPI process to all other MPI processes");

        Bcast_successful = MPI_Bcast(logpress, nrho_ntemp_nye, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast array 'logpress' from the master MPI process to all other MPI processes");

        Bcast_successful = MPI_Bcast(logenergy, nrho_ntemp_nye, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast array 'logenergy' from the master MPI process to all other MPI processes");

        Bcast_successful = MPI_Bcast(dpdrhoe, nrho_ntemp_nye, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast array 'dpdrhoe' from the master MPI process to all other MPI processes");

        Bcast_successful = MPI_Bcast(dpderho, nrho_ntemp_nye, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast array 'dpderho' from the master MPI process to all other MPI processes");

        Bcast_successful = MPI_Bcast(&energy_shift, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (Bcast_successful != MPI_SUCCESS)
            CCTK_ERROR("Could not broadcast variable 'energy_shift' from the master MPI process to all other MPI processes");


        // Synchronize all MPI processes // FIXME: is this really needed?
        MPI_Barrier(MPI_COMM_WORLD);

    #else
        read_arrays_from_EOStable(file);
    #endif



    /* Check whether all logrho, logtemp, ye are sorted into ascending order. If
       they are not, declare failure: in that case, indeed, both binary search
       and temperature root-finding are going to fail.
       FIXME: maybe I can just sort the arrays into ascending order if they are
              not instead of declaring failure? However, all the 3D quantities
              would need to be fixed then...                                    */
    for (i = 1; i < nrho; ++i)
        if (logrho[i] <= logrho[i - 1])
            CCTK_ERROR("Array 'logrho' is not sorted into ascending order.");

    for (j = 1; j < ntemp; ++j)
        if (logtemp[j] <= logtemp[j - 1])
            CCTK_ERROR("Array 'logtemp' is not sorted into ascending order.");

    for (k = 1; k < nye; ++k)
        if (ye[k] <= ye[k - 1])
            CCTK_ERROR("Array 'ye' is not sorted into ascending order.");



    /* Convert log10 to ln, as exp() is way faster than pow(), and convert from
       EOS table units (usually CGS) to geometrized ("Cactus") units. The
       electron fraction ye needs no special treatment.
       N.B.: ln(10^(log10(x)) = log10(x)·ln(10)                                 */
    // FIXME: should I make the following for loops OMP parallel?
    for (i = 0; i < nrho; ++i)
        logrho[i] = ln10*logrho[i] + log_rho_CGStoGeom;

    /* Temperature is already stored in MeV in the EOS table
       => No conversion needed                                                  */
    for (j = 0; j < ntemp; ++j)
        logtemp[j] = ln10*logtemp[j];

    for (kji = 0; kji < nrho_ntemp_nye; ++kji) {
        logpress[kji]  = ln10*logpress[kji]  + log_press_CGStoGeom;
        logenergy[kji] = ln10*logenergy[kji] + log_eps_CGStoGeom;
        dpdrhoe[kji]  *= press_over_rho_CGStoGeom;
        dpderho[kji]  *= press_over_eps_CGStoGeom;
    }

    energy_shift *= eps_CGStoGeom;


    // Define keys corresponding to 3D quantities (declared in RIT_EOS.hh)
    // Pressure key
    press_key.var3D             = logpress;
    press_key.do_exp_log        = true;
    press_key.do_energy_shift   = false;

    // Specific internal energy key
    eps_key.var3D               = logenergy;
    eps_key.do_exp_log          = true;
    eps_key.do_energy_shift     = true;

    // d(press)/d(rho) key
    dpdrhoe_key.var3D           = dpdrhoe;
    dpdrhoe_key.do_exp_log      = false;
    dpdrhoe_key.do_energy_shift = false;

    // d(press)/d(eps) key
    dpderho_key.var3D           = dpderho;
    dpderho_key.do_exp_log      = false;
    dpderho_key.do_energy_shift = false;



    // ***** DEBUG ONLY *****
    /*
    CCTK_REAL xlogrho = -8.12;
    CCTK_REAL xye     = 0.17;
    CCTK_REAL xeps    = exp(-77.24);
    CCTK_REAL eps_out[1];

    CCTK_REAL xlogtemp = logtemp_from_logrho_ye_var3D(xlogrho, xye,
                                                      xeps, eps_key);
    vector<key> which_vars3D{eps_key};
    vars3D_from_logrho_logtemp_ye(xlogrho, xlogtemp, xye,
                                  which_vars3D, eps_out);

    CCTK_VINFO("xeps = %g, xlogtemp = %g, eps_out = %g",
               xeps, xlogtemp, eps_out[0]);
    */

    // ***** END DEBUG *****

    return;
}





// Routine to free memory allocated for EOS table pointers
extern "C" void DeleteTableVars(CCTK_ARGUMENTS) {
    if (logrho == nullptr)
        CCTK_WARN(1, "Pointer 'logrho'    is nullptr, doing nothing. However, this should NOT happen! CHECK WHAT'S GOING ON!!!");
    else {
        delete [] logrho;
        logrho = nullptr;
        CCTK_INFO("Pointer 'logrho'    deleted and set to nullptr");
    }


    if (logtemp == nullptr)
        CCTK_WARN(1, "Pointer 'logtemp'   is nullptr, doing nothing. However, this should NOT happen! CHECK WHAT'S GOING ON!!!");
    else {
        delete [] logtemp;
        logtemp = nullptr;
        CCTK_INFO("Pointer 'logtemp'   deleted and set to nullptr");
    }


    if (ye == nullptr)
        CCTK_WARN(1, "Pointer 'ye'        is nullptr, doing nothing. However, this should NOT happen! CHECK WHAT'S GOING ON!!!");
    else {
        delete [] ye;
        ye = nullptr;
        CCTK_INFO("Pointer 'ye'        deleted and set to nullptr");
    }


    if (logpress == nullptr)
        CCTK_WARN(1, "Pointer 'logpress'  is nullptr, doing nothing. However, this should NOT happen! CHECK WHAT'S GOING ON!!!");
    else {
        delete [] logpress;
        logpress = nullptr;
        CCTK_INFO("Pointer 'logpress'  deleted and set to nullptr");
    }


    if (logenergy == nullptr)
        CCTK_WARN(1, "Pointer 'logenergy' is nullptr, doing nothing. However, this should NOT happen! CHECK WHAT'S GOING ON!!!");
    else {
        delete [] logenergy;
        logenergy = nullptr;
        CCTK_INFO("Pointer 'logenergy' deleted and set to nullptr");
    }


    if (dpdrhoe == nullptr)
        CCTK_WARN(1, "Pointer 'dpdrhoe'   is nullptr, doing nothing. However, this should NOT happen! CHECK WHAT'S GOING ON!!!");
    else {
        delete [] dpdrhoe;
        dpdrhoe = nullptr;
        CCTK_INFO("Pointer 'dpdrhoe'   deleted and set to nullptr");
    }


    if (dpderho == nullptr)
        CCTK_WARN(1, "Pointer 'dpderho'   is nullptr, doing nothing. However, this should NOT happen! CHECK WHAT'S GOING ON!!!");
    else {
        delete [] dpderho;
        dpderho = nullptr;
        CCTK_INFO("Pointer 'dpderho'   deleted and set to nullptr");
    }


    return;
}
