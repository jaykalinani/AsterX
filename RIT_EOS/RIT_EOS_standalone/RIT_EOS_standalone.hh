#ifndef RIT_EOS_STANDALONE_HH
#define RIT_EOS_STANDALONE_HH

#include <string>

using namespace std;



extern "C" {
    // ############################ I/O ROUTINES ###############################
    // Routine to get the dimensions of the table
    void RIT_EOS_ReadTable(const string     EOS_TableName,
                                         int        &nrho,
                                         int        &ntemp,
                                         int        &nye,
                                         double    *&logrho,
                                         double    *&logtemp,
                                         double    *&ye,
                                         double  ***&logpress,
                                         double  ***&logenergy,
                                         double  ***&entropy);

    /* Routine to free the memory allocated for logrho, logtemp, ye, logpress,
       logenergy, entropy by routine RIT_EOS_ReadTable                  */
    void RIT_EOS_DeleteTableVars(double   *&logrho,
                                         double   *&logtemp,
                                         double   *&ye,
                                         double ***&logpress,
                                         double ***&logenergy,
                                         double ***&entropy);





    // ############### ROUTINES TO GET INTERPOLATION BOUNDS ####################
    /* Routines to get the tabulated values of a some quantity which are closest
       to some chosen value for that quantity, which must be in the table's
       range (otherwise execution is aborted)                                   */
    void RIT_EOS_GetBounds_1Dvars(const double  &x,       // xlogrho, xlogtemp, xye
                                                double *&_1Darr,  // logrho,  logtemp,  ye   // FIXME: 'const' not accepted
                                          const int     &n,       // nrho,    ntemp,    nye
                                                int    (&bounds)[2]);

    void RIT_EOS_GetBounds_3Dvars_logrho_ye_tabulated(const double    &x,
                                                                    double ***&_3Darr,
                                                              const int       &index_logrho,
                                                              const int       &index_ye,
                                                              const int       &nrho,
                                                              const int       &ntemp,
                                                              const int       &nye,
                                                                    int      (&bounds)[2]);





    // ############### ROUTINES TO INTERPOLATE THE EOS TABLE ###################

    /* Routine to retrieve a value from a 3D array for values of logrho,
       logtemp, ye which are not tabulated; the cheapest way to do this is
       trilinear interpolation.                                                 */
    double RIT_EOS_interp_3Darr_from_logrho_logtemp_ye(      double ***&_3Darr,
                                                               const double    &xlogrho,
                                                               const double    &xlogtemp,
                                                               const double    &xye,
                                                                     double   *&logrho,
                                                                     double   *&logtemp,
                                                                     double   *&ye,
                                                               const int       &nrho,
                                                               const int       &ntemp,
                                                               const int       &nye);//,
                                                           /*  const int      (&bounds_logrho)[2],
                                                               const int      (&bounds_logtemp)[2],
                                                               const int      (&bounds_ye)[2]);*/

    // Routine to retrieve logtemp knowing logrho, ye, and one 3D array
    double RIT_EOS_get_logtemp_from_logrho_ye_3Dvar(const double    &xlogrho,
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
                                                            const int       &itmax_rootfinding);
}

#endif  // RIT_EOS_STANDALONE_HH
