#include <iostream>
#include <string>
#include <cstdlib>  // For exit()
#include <cmath>    // For fabs()
#include <hdf5/serial/hdf5.h>
#include "RIT_EOS_standalone.hh"

using namespace std;





// ############################ I/O ROUTINES ###################################

/* Routine to get the dimensions of the table

   ***** WARNING *****
   Remember to delete the arrays allocated by this routine (the simplest way to
   do this is to call the routine RIT_EOS_DeleteTableVars)

   ***** IMPORTANT NOTES ****
   This function accepts *REFERENCES TO POINTERS* as input. The reasons for this
   choice are the following:

   1. passing by reference (C++ only!) is more efficient than passing by value,
      because no copy is made of the object passed to the function;
   2. a reference is an alias of the the variable in the caller function, not a
      copy of it, so modifying the referenced object in the called function has
      an effect on the object in the caller function (this would not be the case
      when passing by value).

   A few remarks:
   1. it would be 'cleaner' to pass a reference to an ARRAY instead of a
      reference to a pointer (something like "double (&logrho)[nrho]", with nrho
      known at COMPILE time). However, EOS tables may be very large and arrays
      are allocated on the STACK; but allocating large amounts of memory on the
      stack is unadvisable and not portable (can lead to stack overflows on some
      machines). A pointer is instead allocated (dynamically) on the HEAP and
      there is no threat of memory issues (provided that the available RAM is
      not exceeded). On the other hand, heap allocation is slower than stack
      allocation, but it is definitely worth it here;
   2. a simpler strategy would be to have GLOBAL POINTERS, but this is poor and
      bad programming technique. Non-const global variables can be accessed and
      modified from any part of the code and are very unsafe;
   3. it would be good practice to use VECTORS instead of pointers, but the
      following routine is meant to be a scheduled function in the Einstein
      Toolkit version of RIT_EOS and scheduled functions do not work
      with vectors.                                                             */

extern "C" void RIT_EOS_ReadTable(const string     EOS_TableName,
                                                int       &nrho,
                                                int       &ntemp,
                                                int       &nye,
                                                double   *&logrho,
                                                double   *&logtemp,
                                                double   *&ye,
                                                double ***&logpress,
                                                double ***&logenergy,
                                                double ***&entropy) {
    // Declaration of variables
    int     j, k;
    int     k_ntemp, j_nrho, k_ntemp_nrho;
    herr_t  status;
    hid_t   file_id;
    hid_t   dataset_id;

    // Open HDF5 EOS table file
    file_id = H5Fopen(EOS_TableName.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

    // Handle errors in opening EOS table file
    if (file_id < 0) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_CheckTableDims\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Unable to open EOS table " << EOS_TableName
             << endl << endl;
        exit(1);
    }


    // Get the number of points in logrho, logtemp, ye from the table
    dataset_id = H5Dopen2(file_id, "pointsrho", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nrho);
    status     = H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "pointstemp", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &ntemp);
    status     = H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "pointsye", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nye);
    status     = H5Dclose(dataset_id);



    // Check if the dimensions of the EOS table make sense (they should)
    // FIXME: change 'cout' + 'exit' to CCTK_WARNING or similar
    if (nrho < 1) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
               "\nRoutine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
               "\nLine      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
               "\n-> Problem with the EOS table: number of points in array 'logrho' is lower than 1."
             << endl << endl;
        exit(1);
    }

    if (ntemp < 1) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
               "\nRoutine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
               "\nLine      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
               "\n-> Problem with the EOS table: number of points in array 'logtemp' is lower than 1."
             << endl << endl;
        exit(1);
    }

    if (nye < 1) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
               "\nRoutine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
               "\nLine      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
               "\n-> Problem with the EOS table: number of points in array 'ye' is lower than 1."
             << endl << endl;
        exit(1);
    }



    // Allocate memory for the quantities to be read from the table
    logrho          = new double[nrho];
    logtemp         = new double[ntemp];
    ye              = new double[nye];

    logpress        = new double**[nye];
    logpress[0]     = new double*[nye*ntemp];
    logpress[0][0]  = new double[nye*ntemp*nrho];

    logenergy       = new double**[nye];
    logenergy[0]    = new double*[nye*ntemp];
    logenergy[0][0] = new double[nye*ntemp*nrho];

    entropy         = new double**[nye];
    entropy[0]      = new double*[nye*ntemp];
    entropy[0][0]   = new double[nye*ntemp*nrho];

    for (k = 0; k < nye; ++k) {
        k_ntemp      = k*ntemp;
        logpress[k]  = logpress[0]  + k_ntemp;
        logenergy[k] = logenergy[0] + k_ntemp;
        entropy[k]   = entropy[0]   + k_ntemp;
        for (j = 0; j < ntemp; ++j) {
            j_nrho          = j*nrho;
            k_ntemp_nrho    = k*ntemp*nrho;
            logpress[k][j]  = logpress[0][0]  + (k_ntemp_nrho) + j_nrho;
            logenergy[k][j] = logenergy[0][0] + (k_ntemp_nrho) + j_nrho;
            entropy[k][j]   = entropy[0][0]   + (k_ntemp_nrho) + j_nrho;
        }
    }


    // Check if pointers have been allocated correctly
    // FIXME: change 'printf' + 'exit' to CCTK_WARNING(0, ...) or similar
    if (logrho == nullptr) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Failed to allocate memory for array 'logrho'"
             << endl << endl;
        exit(1);
    }

    if (logtemp == nullptr) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Failed to allocate memory for array 'logtemp'"
             << endl << endl;
        exit(1);
    }

    if (ye == nullptr) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Failed to allocate memory for array 'ye'"
             << endl << endl;
        exit(1);
    }

    if (logpress == nullptr) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Failed to allocate memory for array 'logpress'"
             << endl << endl;
        exit(1);
    }

    if (logenergy == nullptr) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Failed to allocate memory for array 'logenergy'"
             << endl << endl;
        exit(1);
    }

    if (entropy == nullptr) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_ReadTable\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Failed to allocate memory for array 'entropy'"
             << endl << endl;
        exit(1);
    }



    // Read datasets from the table
    dataset_id = H5Dopen2(file_id, "logrho", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, logrho);
    status     = H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "logtemp", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, logtemp);
    status     = H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "ye", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ye);
    status     = H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "logpress", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, logpress[0][0]);
    status     = H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "logenergy", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, logenergy[0][0]);
    status     = H5Dclose(dataset_id);

    dataset_id = H5Dopen2(file_id, "entropy", H5P_DEFAULT);
    status     = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, entropy[0][0]);
    status     = H5Dclose(dataset_id);

    status     = H5Fclose(file_id);
    return;
}





/* Routine to free the memory allocated for logrho, logtemp, ye, logpress,
   logenergy, entropy by routine RIT_EOS_ReadTable                      */
extern "C" void RIT_EOS_DeleteTableVars(double   *&logrho,
                                                double   *&logtemp,
                                                double   *&ye,
                                                double ***&logpress,
                                                double ***&logenergy,
                                                double ***&entropy) {
    // FIXME: change 'cout' to CCTK_INFO or similar
    if (logrho == nullptr)
        cout << endl << "\033[0;33m\033[1mWARNING\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_DeleteTableVars\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'logrho' either not initialised or already freed. Doing nothing."
             << endl << endl;
    else {
        delete [] logrho;
        logrho = nullptr;
        cout << endl << "\033[0;32m\033[1mINFO\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'logrho' freed and set to 'nullptr'" << endl << endl;
    }

    if (logtemp == nullptr)
        cout << endl << "\033[0;33m\033[1mWARNING\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_DeleteTableVars\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'logtemp' either not initialised or already freed. Doing nothing."
             << endl << endl;
    else {
        delete [] logtemp;
        logtemp = nullptr;
        cout << endl << "\033[0;32m\033[1mINFO\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'logtemp' freed and set to 'nullptr'" << endl << endl;
    }

    if (ye == nullptr)
        cout << endl << "\033[0;33m\033[1mWARNING\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_DeleteTableVars\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'ye' either not initialised or already freed. Doing nothing."
             << endl << endl;
    else {
        delete [] ye;
        ye = nullptr;
        cout << endl << "\033[0;32m\033[1mINFO\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'ye' freed and set to 'nullptr'" << endl << endl;
    }

    if (logpress == nullptr)
        cout << endl << "\033[0;33m\033[1mWARNING\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_DeleteTableVars\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'logpress' either not initialised or already freed. Doing nothing."
             << endl << endl;
    else {
        delete [] logpress;
        logpress = nullptr;
        cout << endl << "\033[0;32m\033[1mINFO\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'logpress' freed and set to 'nullptr'" << endl << endl;
    }

    if (logenergy == nullptr)
        cout << endl << "\033[0;33m\033[1mWARNING\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_DeleteTableVars\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'logenergy' either not initialised or already freed. Doing nothing."
             << endl << endl;
    else {
        delete [] logenergy;
        logenergy = nullptr;
        cout << endl << "\033[0;32m\033[1mINFO\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'logenergy' freed and set to 'nullptr'" << endl << endl;
    }

    if (entropy == nullptr)
        cout << endl << "\033[0;33m\033[1mWARNING\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_DeleteTableVars\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'entropy' either not initialised or already freed. Doing nothing."
             << endl << endl;
    else {
        delete [] entropy;
        entropy = nullptr;
        cout << endl << "\033[0;32m\033[1mINFO\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Array 'entropy' freed and set to 'nullptr'" << endl << endl;
    }


    return;
}

// #############################################################################










// ############### ROUTINES TO GET INTERPOLATION BOUNDS ########################

/* Routines to get the tabulated values of a some quantity which are closest to
   some chosen value for that quantity, which must be in the table's range
   (otherwise execution is aborted)                                             */
extern "C" void RIT_EOS_GetBounds_1Dvar(const double &x,             // xlogrho, xlogtemp, xye
                                                      double      *&_1Darr,  // logrho,  logtemp,  ye   // FIXME: 'const' not accepted
                                				      const int    &n,       // nrho,    ntemp,    nye
                                				      int         (&bounds)[2]) {
    // Preliminary check
    if (x < _1Darr[0] || x > _1Darr[n - 1]) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Desired point (" << x << ") is out of bounds."
             << endl << endl;
        exit(1);
	}

   /* FIXME: 1) check if '_1Darr' is a null pointer? But this can be slow if this
                routine is called many times;
             2) checks on the value of n?                                       */


    // Declaration of variables
    int ilow  = 0;
    int iup   = n - 1;
    int ihalf;

    // Locate x inside arr using binary search
    while (ilow < iup - 1) {
        ihalf = (ilow + iup)/2;  // !!! INTEGER division !!!
        if (x >= _1Darr[ihalf])  ilow = ihalf;
        else                    iup  = ihalf;
    }

    // Set up bounds for x
    bounds[0] = ilow;
    bounds[1] = iup;

    return;
}





extern "C" void RIT_EOS_GetBounds_3Dvar_logrho_ye_tabulated(const double    &x,
                                                                          double ***&_3Darr,
                                                                    const int       &index_logrho,
                                                                    const int       &index_ye,
                                                                    const int       &nrho,
                                                                    const int       &ntemp,
                                                                    const int       &nye,
                                                                          int      (&bounds)[2]) {
    // Preliminary checks
    if (index_logrho < 0 || index_logrho > nrho - 1) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds_3Dvar_logrho_ye_tabulated\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Desired interpolation point is out of bounds in 'logrho' direction (index_logrho = "
                     << index_logrho << ", index_logrho_min = 0, index_logrho_max = " << nrho - 1
             << endl << endl;
        exit(1);
	}

    if (index_ye < 0 || index_ye > nye - 1) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds_3Dvar_logrho_ye_tabulated\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Desired interpolation point is out of bounds in 'ye' direction (index_ye = "
                     << index_ye << ", index_ye_min = 0, index_ye_max = " << nye - 1
             << endl << endl;
        exit(1);
	}

    if (x < _3Darr[index_ye][0][index_logrho] ||
        x > _3Darr[index_ye][ntemp - 1][index_logrho]) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_GetBounds_3Dvar_logrho_ye_tabulated\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "Desired interpolation point is out of bounds along the 3D variable in 'logtemp' direction (x3Dvar = "
                     << _3Darr[index_ye][0][index_logrho]
                     << ", _3Darr[" << index_ye << "][0][" << index_logrho << "] = "
                     << _3Darr[index_ye][0][index_logrho] << ","
                     << "_3Darr[index_ye][ntemp - 1][index_logrho] = "
                     << _3Darr[index_ye][ntemp - 1][index_logrho] << ")."
             << endl << endl;
        exit(1);
    }

   /* FIXME: 1) check if '_3Darr' is a null pointer? But this can be slow if this
                routine is called many times;
             2) checks on the values of nrho, ntemp, nye?                       */


    // Declaration of variables
    int jlow  = 0;
    int jup   = ntemp - 1;
    int jhalf;

    // Locate x inside array arr[index_ye][:][index_logrho] using binary search
    while (jlow < jup - 1) {
        jhalf = (jlow + jup)/2;  // !!! INTEGER division !!!
        if (x >= _3Darr[index_ye][jhalf][index_logrho])  jlow = jhalf;
        else                                             jup  = jhalf;
    }

    // Set up bounds for x
    bounds[0] = jlow;
    bounds[1] = jup;

    return;
}

// #############################################################################










// ############### ROUTINES TO INTERPOLATE THE EOS TABLE #######################

/* Routine to retrieve a value from a 3D array for values of logrho, logtemp, ye
   which are not tabulated; the cheapest way to do this is trilinear
   interpolation.                                                               */
extern "C" double RIT_EOS_interp_3Darr_from_logrho_logtemp_ye(      double ***&_3Darr,
                                                                      const double    &xlogrho,
                                                                      const double    &xlogtemp,
                                                                      const double    &xye,
                                                                            double   *&logrho,
                                                                            double   *&logtemp,
                                                                            double   *&ye,
                                                                      const int       &nrho,
                                                                      const int       &ntemp,
                                                                      const int       &nye) {
    // Preliminary checks
    if (xlogrho < logrho[0] || xlogrho > logrho[nrho - 1]) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_interp_3Dvars_from_logrho_logtemp_ye\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Desired interpolation point is out of bounds in 'logrho' direction (xlogrho = "
                     << xlogrho << ", logrho[0] = "
                     << logrho[0] << ", logrho[" << nrho - 1 << "] = "
                     << logrho[nrho - 1] << ")."
             << endl << endl;
        exit(1);
	}

    if (xlogtemp < logtemp[0] || xlogtemp > logtemp[ntemp - 1]) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_interp_3Dvars_from_logrho_logtemp_ye\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Desired interpolation point is out of bounds in 'logtemp' direction (xlogtemp = "
                     << xlogtemp << ", logtemp[0] = "
                     << logtemp[0] << ", logtemp[" << ntemp - 1 << "] = "
                     << logtemp[ntemp - 1] << ")."
             << endl << endl;
        exit(1);
	}

    if (xye < ye[0] || xye > ye[nye - 1]) {
        cout << endl << "\033[0;31m\033[1mERROR\033[0m"
             << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
             << endl << "Routine   \033[0;36m\033[1mRIT_EOS_interp_3Dvars_from_logrho_logtemp_ye\033[0m"
             << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
             << endl << "-> Desired interpolation point is out of bounds in 'ye' direction (xye = "
                     << xye << ", ye[0] = "
                     << ye[0] << ", ye[" << nye - 1 << "] = "
                     << ye[nye - 1] << ")."
             << endl << endl;
        exit(1);
	}


    // Declaration of variables
    int     bounds_logrho[2];
    int     bounds_logtemp[2];
    int     bounds_ye[2];
    int     ilow, iup, jlow, jup, klow, kup;
    double  deltalow, deltaup, deltatot_inv;
    double  aux1, aux2, aux3, aux4, aux5, aux6;

    // Get bounds for logrho, logtemp and ye
    RIT_EOS_GetBounds_1Dvar(xlogrho,  logrho,  nrho,  bounds_logrho);
    RIT_EOS_GetBounds_1Dvar(xlogtemp, logtemp, ntemp, bounds_logtemp);
    RIT_EOS_GetBounds_1Dvar(xye,      ye,      nye,   bounds_ye);

    // Set up auxiliary variables for interpolation
    ilow = bounds_logrho[0];
    iup  = bounds_logrho[1];
    jlow = bounds_logtemp[0];
    jup  = bounds_logtemp[1];
    klow = bounds_ye[0];
    kup  = bounds_ye[1];

    // Partial interpolations for fixed logrho
    deltalow     = xlogrho - logrho[ilow];
    deltaup      = xlogrho - logrho[iup];
    deltatot_inv = 1./(1.*(logrho[iup] - logrho[ilow]));
    aux1         = (_3Darr[klow][jlow][iup]*deltalow - _3Darr[klow][jlow][ilow]*deltaup)*deltatot_inv;
    aux2         = (_3Darr[klow][jup][iup] *deltalow - _3Darr[klow][jup][ilow] *deltaup)*deltatot_inv;
    aux3         = (_3Darr[kup][jlow][iup] *deltalow - _3Darr[kup][jlow][ilow] *deltaup)*deltatot_inv;
    aux4         = (_3Darr[kup][jup][iup]  *deltalow - _3Darr[kup][jup][ilow]  *deltaup)*deltatot_inv;

    // Partial interpolations for fixed logtemp
    deltalow     = xlogtemp - logtemp[jlow];
    deltaup      = xlogtemp - logtemp[jup];
    deltatot_inv = 1./(1.*(logtemp[jup] - logtemp[jlow]));
    aux5         = (aux2*deltalow - aux1*deltaup)*deltatot_inv;
    aux6         = (aux4*deltalow - aux3*deltaup)*deltatot_inv;

    // Final interpolation (fixed ye)
    return (aux6*(xye - ye[klow]) - aux5*(xye - ye[kup]))/(1.*(ye[kup] - ye[klow]));
}










// Routine to retrieve logtemp knowing logrho, ye, and one 3D array
extern "C" double RIT_EOS_get_logtemp_from_logrho_ye_3Dvar(const double    &xlogrho,
                                                                   const double    &xye,
                                                                   const double    &x3Dvar,
                                                                         double ***&_3Darr,
                                                                         double   *&logrho,
                                                                         double   *&logtemp,
                                                                         double   *&ye,
                                                                   const int       &nrho,
                                                                   const int       &ntemp,
                                                                   const int       &nye,
                                                                   const double    &eps_rootfinding,
                                                                   const int       &itmax_rootfinding) {
    // Declaration of variables
    int    bounds_logrho[2];
    int    bounds_ye[2];
    int    bounds_3Darr_tabulated_logrho_ye_1[2];
    int    bounds_3Darr_tabulated_logrho_ye_2[2];
    int    bounds_3Darr_tabulated_logrho_ye_3[2];
    int    bounds_3Darr_tabulated_logrho_ye_4[2];
    double xlogtemp, xlogtemp_bound_low, xlogtemp_bound_up;
    double f, flow, fup;
    double _3Dvar_jlow;
    int    j, jlow, jup;
    int    is_3Dvar_jlow_LT_x3Dvar;
    int    n;


    // Get bounds for logrho and ye
    RIT_EOS_GetBounds_1Dvar(xlogrho, logrho, nrho, bounds_logrho);
    RIT_EOS_GetBounds_1Dvar(xye,     ye,     nye,  bounds_ye);

    /* Get bounds for _3Darr on four different lines of varying
       logtemp:
       1. for fixed   logrho(bounds_logrho[0]),
                      ye(bounds_ye[0]);
       2. for fixed   logrho(bounds_logrho[0]),
                      ye(bounds_ye[1]);
       3. for fixed   logrho(bounds_logrho[1]),
                      ye(bounds_ye[0]);
       4. for fixed   logrho(bounds_logrho[1]),
                      ye(bounds_ye[1]).                                         */

    RIT_EOS_GetBounds_3Dvar_logrho_ye_tabulated(x3Dvar, _3Darr,
        bounds_logrho[0], bounds_ye[0], nrho, ntemp, nye,
        bounds_3Darr_tabulated_logrho_ye_1);

    RIT_EOS_GetBounds_3Dvar_logrho_ye_tabulated(x3Dvar, _3Darr,
        bounds_logrho[0], bounds_ye[1], nrho, ntemp, nye,
        bounds_3Darr_tabulated_logrho_ye_2);

    RIT_EOS_GetBounds_3Dvar_logrho_ye_tabulated(x3Dvar, _3Darr,
        bounds_logrho[1], bounds_ye[0], nrho, ntemp, nye,
        bounds_3Darr_tabulated_logrho_ye_3);

    RIT_EOS_GetBounds_3Dvar_logrho_ye_tabulated(x3Dvar, _3Darr,
        bounds_logrho[1], bounds_ye[1], nrho, ntemp, nye,
        bounds_3Darr_tabulated_logrho_ye_4);


    // Set bounds to bracket xlogtemp
    jlow = min(min(bounds_3Darr_tabulated_logrho_ye_1[0],
                   bounds_3Darr_tabulated_logrho_ye_2[0]),
               min(bounds_3Darr_tabulated_logrho_ye_3[0],
                   bounds_3Darr_tabulated_logrho_ye_4[0]));
    jup  = max(max(bounds_3Darr_tabulated_logrho_ye_1[1],
                   bounds_3Darr_tabulated_logrho_ye_2[1]),
               max(bounds_3Darr_tabulated_logrho_ye_3[1],
                   bounds_3Darr_tabulated_logrho_ye_4[1]));

    if (jup == jlow + 1) {
        xlogtemp_bound_low = logtemp[jlow];
        xlogtemp_bound_up  = logtemp[jup];
    }

    else {
        j                       = jlow + 1;
        _3Dvar_jlow             = RIT_EOS_interp_3Darr_from_logrho_logtemp_ye(
                                      _3Darr, xlogrho, logtemp[jlow], xye,
                                      logrho, logtemp, ye, nrho, ntemp, nye);
        is_3Dvar_jlow_LT_x3Dvar = (_3Dvar_jlow < x3Dvar); // 1 if _3Dvar_jlow < x3Dvar, 0 otherwise

        while ((RIT_EOS_interp_3Darr_from_logrho_logtemp_ye(
                    _3Darr, xlogrho, logtemp[j], xye,
                    logrho, logtemp, ye, nrho, ntemp, nye) < x3Dvar)
                == is_3Dvar_jlow_LT_x3Dvar)
            ++j;

        xlogtemp_bound_low = logtemp[j - 1];
        xlogtemp_bound_up  = logtemp[j];
    }


    /* Find xlogtemp such that _3Dvar(xye, xlogtemp, xlogrho) = x3Dvar
       using bisection                                                          */
    for (n = 0; n < itmax_rootfinding; ++n) {
        xlogtemp  = 0.5*(xlogtemp_bound_up + xlogtemp_bound_low);
        if (fabs(xlogtemp_bound_up - xlogtemp_bound_low) < eps_rootfinding)
            return xlogtemp;

        f = RIT_EOS_interp_3Darr_from_logrho_logtemp_ye(
                _3Darr, xlogrho, xlogtemp, xye,
                logrho, logtemp, ye, nrho, ntemp, nye);
        if      (f < x3Dvar)  xlogtemp_bound_low = xlogtemp;
        else if (f > x3Dvar)  xlogtemp_bound_up  = xlogtemp;
        else    return xlogtemp;
    }


    cout << endl << "\033[0;31m\033[1mERROR\033[0m"
         << endl << "From file \033[0;36m\033[1m" << __FILE__ << "\033[0m"
         << endl << "Routine   \033[0;36m\033[1mRIT_EOS_get_logtemp_from_logrho_ye_3Dvar\033[0m"
         << endl << "Line      \033[0;36m\033[1m" << __LINE__ << "\033[0m"
         << endl << "-> Bisection failed (more than " << itmax_rootfinding
                 << " iterations done without converging to precision "
                 << eps_rootfinding << "."
         << endl << endl;
    exit(1);
}



