#include <iostream>
#include <string>
#include "RIT_EOS_standalone.hh"


// ######################### USER-DEFINED PARAMETERS ###########################

#define XLOGRHO          9.1122
#define XLOGTEMP         0.8832
#define XYE              0.1791
#define XLOGPRESS        32.22
#define XLOGENERGY       19.77
#define XENTROPY         46.653
#define EPS_ROOTFINDING  1.e-10
#define ITMAX_BISECTION  200
#define EOS_TABLENAME    "LS220star_0000_rho391_temp163_ye66.h5"

// #############################################################################





int main() {
    // Declaration of variables
    const string EOS_TableName = EOS_TABLENAME;
    int          nrho,        ntemp,        nye;
    double      *logrho,     *logtemp,     *ye;
    double    ***logpress, ***logenergy, ***entropy;
/*    int          bounds_logrho[2];
    int          bounds_logtemp[2];
    int          bounds_ye[2];*/
    int          bounds_logpress_tabulated_logrho_ye[2];
    int          bounds_logenergy_tabulated_logrho_ye[2];
    int          bounds_entropy_tabulated_logrho_ye[2];
    double       retrieved_logtemp;


    // Read the EOS table
    RIT_EOS_ReadTable(EOS_TableName, nrho, ntemp, nye, logrho, logtemp,
                              ye, logpress, logenergy, entropy);


    // Get bounds for XLOGRHO, XLOGTEMP, XYE
/*    RIT_EOS_GetBounds_1Dvars(XLOGRHO,  logrho,  nrho,  bounds_logrho);
    RIT_EOS_GetBounds_1Dvars(XLOGTEMP, logtemp, ntemp, bounds_logtemp);
    RIT_EOS_GetBounds_1Dvars(XYE,      ye,      nye,   bounds_ye);*/


    /* Get bounds for XLOGPRESS, XLOGENERGY, XENTROPY at fixed
       logrho and ye                                                            */
    /*RIT_EOS_GetBounds_3Dvars_logrho_ye_tabulated(XLOGPRESS, logpress,
        bounds_logrho[0], bounds_ye[0], nrho, ntemp, nye,
        bounds_logpress_tabulated_logrho_ye);

    RIT_EOS_GetBounds_3Dvars_logrho_ye_tabulated(XLOGENERGY, logenergy,
        bounds_logrho[0], bounds_ye[0], nrho, ntemp, nye,
        bounds_logenergy_tabulated_logrho_ye);

    RIT_EOS_GetBounds_3Dvars_logrho_ye_tabulated(XENTROPY, entropy,
        bounds_logrho[0], bounds_ye[0], nrho, ntemp, nye,
        bounds_entropy_tabulated_logrho_ye);*/





    // ######################### RETRIEVE TEMPERATURE ##########################
    // Retrieve logtemp from logrho, ye, logpress
    retrieved_logtemp = RIT_EOS_get_logtemp_from_logrho_ye_3Dvar(
                            XLOGRHO, XYE, XLOGPRESS,
                            logpress, logrho, logtemp, ye,
                            nrho, ntemp, nye,
                            EPS_ROOTFINDING, ITMAX_BISECTION);

    cout << endl << "logpress(" << XYE << ", logtemp, " << XLOGRHO << ") = "
                 << XLOGPRESS << " for logtemp = " << retrieved_logtemp << endl;

//    RIT_EOS_GetBounds_1Dvars(retrieved_logtemp, logtemp, ntemp, bounds_logtemp);

    cout << "CHECK: logpress(" << XYE << ", " << retrieved_logtemp << ", " << XLOGRHO << ") = "
         << RIT_EOS_interp_3Darr_from_logrho_logtemp_ye(
                logpress, XLOGRHO, retrieved_logtemp, XYE, logrho, logtemp, ye,
                nrho, ntemp, nye)//, bounds_logrho, bounds_logtemp, bounds_ye)
         << endl << endl;


    // Retrieve logtemp from logrho, ye, logenergy
    retrieved_logtemp = RIT_EOS_get_logtemp_from_logrho_ye_3Dvar(
                            XLOGRHO, XYE, XLOGENERGY,
                            logenergy, logrho, logtemp, ye,
                            nrho, ntemp, nye,
                            EPS_ROOTFINDING, ITMAX_BISECTION);

    cout << endl << "logenergy(" << XYE << ", logtemp, " << XLOGRHO << ") = "
                 << XLOGENERGY << " for logtemp = " << retrieved_logtemp << endl;

//    RIT_EOS_GetBounds_1Dvars(retrieved_logtemp, logtemp, ntemp, bounds_logtemp);

    cout << "CHECK: logenergy(" << XYE << ", " << retrieved_logtemp << ", " << XLOGRHO << ") = "
         << RIT_EOS_interp_3Darr_from_logrho_logtemp_ye(
                logenergy, XLOGRHO, retrieved_logtemp, XYE, logrho, logtemp, ye,
                nrho, ntemp, nye)//, bounds_logrho, bounds_logtemp, bounds_ye)
         << endl << endl;


    // Retrieve logtemp from logrho, ye, entropy
    retrieved_logtemp = RIT_EOS_get_logtemp_from_logrho_ye_3Dvar(
                        XLOGRHO, XYE, XENTROPY,
                        entropy, logrho, logtemp, ye,
                        nrho, ntemp, nye,
                        EPS_ROOTFINDING, ITMAX_BISECTION);

    cout << endl << "entropy(" << XYE << ", logtemp, " << XLOGRHO << ") = "
                 << XENTROPY << " for logtemp = " << retrieved_logtemp << endl;

//    RIT_EOS_GetBounds_1Dvars(retrieved_logtemp, logtemp, ntemp, bounds_logtemp);

    cout << "CHECK: entropy(" << XYE << ", " << retrieved_logtemp << ", " << XLOGRHO << ") = "
         << RIT_EOS_interp_3Darr_from_logrho_logtemp_ye(
                entropy, XLOGRHO, retrieved_logtemp, XYE, logrho, logtemp, ye,
                nrho, ntemp, nye)//, bounds_logrho, bounds_logtemp, bounds_ye)
         << endl << endl;

    // #########################################################################



    // Free the memory allocated by RIT_EOS_ReadTable
    RIT_EOS_DeleteTableVars(logrho, logtemp, ye, logpress, logenergy, entropy);

    return 0;
}
