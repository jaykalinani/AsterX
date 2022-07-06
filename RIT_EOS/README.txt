Cactus Code Thorn RIT_EOS
Author(s)    : Lorenzo Ennoggi <lorenzo.ennoggi@gmail.com> <le8016@rit.edu>
Maintainer(s): Lorenzo Ennoggi <lorenzo.ennoggi@gmail.com> <le8016@rit.edu>
Licence      : GPL
--------------------------------------------------------------------------


1. Purpose of this thorn

This thorn provides routines to retrieve a number of thermodynamic quantities
from an EOS table in Stellarcollapse format (see
<https://stellarcollapse.org/microphysics.html> for an overview and the actual
tables and refer to the documentation --
<https://bitbucket.org/andschn/sroeos/src/master/User_Guide/User_Guide.pdf> --
for futher details).



2. How to use this thorn

A few parameters must be set in your parfile:
    RIT_EOS::EOStable_path   = <Absolute/path/to/EOS/table>
    RIT_EOS::EOS_type        = 1 #2    # 1 for polytropic EOS, 2 for tabulated EOS
    RIT_EOS::bisection_eps   = 1.e-10  # Trade-off between accuracy and speed
    RIT_EOS::bisection_maxit = 200     # Should be a reasonable value
    RIT_EOS::poly_K          = 100.    # Polytropic constant K in c=G=1 units
    RIT_EOS::poly_Gamma      = 2.      # Polytropic index Gamma in c=G=1 units
