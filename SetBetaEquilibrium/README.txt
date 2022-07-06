Cactus Code Thorn SetBetaEquilibrium
Author(s)    : Lorenzo Ennoggi <lorenzo.ennoggi@gmail.com> <le8016@rit.edu>
Maintainer(s): Lorenzo Ennoggi <lorenzo.ennoggi@gmail.com> <le8016@rit.edu>
Licence      : GPL
--------------------------------------------------------------------------



1. Purpose of this thorn

This thorn reads in beta-equilibrated temperature and electron fraction from
initial data in Compose or Lorene format.

********************************************************************************
N.B.: make sure you are using the same EOS which has been used to generate the
      initial data!
********************************************************************************



2. How to use this thorn

A few parameters must be set in your parfile:

    SetBetaEquilibrium::ID_format = "Compose" or "Lorene"
    SetBetaEquilibrium::ID_path   = </Path/to/EOS/slice>                                            in Lorene format or
                                    </Path/to/EOS/slice>/<Name_of_the_EOS_files_without_extension>  in Compose format

    SetBetaEquilibrium::Lorene_slice_filename = <Filename of the initial data slice (Lorene format only)>
    SetBetaEquilibrium::Ye_filename           = <Filename of the initial electron fraction table (Lorene format only)>

    SetBetaEquilibrium::get_EOS_eps_from_beta_eq_ID = "yes" #"no"
    SetBetaEquilibrium::SetPressure                 = "yes" #"no"

    EOS_Omni::nuceos_table_name = </path/to/EOS/table> (in StellarCollapse format)

    HydroBase::initial_temperature = "SetBetaEquilibrium"
    HydroBase::initial_Y_e         = "SetBetaEquilibrium"
    HydroBase::initial_entropy     = "SetBetaEquilibrium"
