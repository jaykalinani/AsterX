#include <string>
#include <fstream>

#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

/*  ***** NOTICE *****
    Including <filesystem> requires std=c++17 (or std=gnu++17) as a compilation
    option (i.e., set CXX_FLAGS = std=c++17 or CXX_FLAGS = std=gnu++17).        */
//#include <experimental/filesystem>  // FIXME: can't link on Frontera

using std::ifstream;
using std::string;



extern "C" void CheckParameters (CCTK_ARGUMENTS) {
    DECLARE_CCTK_ARGUMENTS_CHECKED(CheckParameters);
    DECLARE_CCTK_PARAMETERS;

    // Declaration of variables
    string   ID_path_string = ID_path;
    string   filepath;
    ifstream file;


    // Basic consistency check 1
    if (CCTK_EQUALS(initial_temperature, "SetBetaEquilibrium")) {
        if (not(CCTK_EQUALS(initial_Y_e, "SetBetaEquilibrium")))
            CCTK_PARAMWARN("Since HydroBase::initial_temperature = \"SetBetaEquilibrium\", please also set HydroBase::initial_Y_e = \"SetBetaEquilibrium\".");

        if (CCTK_EQUALS(ID_path, ""))
            CCTK_PARAMWARN("Since HydroBase::initial_temperature = \"SetBetaEquilibrium\", please also set SetBetaEquilibrium::ID_path = </path/to/EOS/slice> (Lorene format) or SetBetaEquilibrium::ID_path = </path/to/EOS/slice>/eos (Compose format).");
    }

    // Basic consistency check 2
    if (CCTK_EQUALS(initial_Y_e, "SetBetaEquilibrium")) {
        if (not(CCTK_EQUALS(initial_temperature, "SetBetaEquilibrium")))
            CCTK_PARAMWARN("Since HydroBase::initial_Y_e = \"SetBetaEquilibrium\", please also set HydroBase::initial_temperature = \"SetBetaEquilibrium\".");

        if (CCTK_EQUALS(ID_path, ""))
            CCTK_PARAMWARN("Since HydroBase::initial_Y_e = \"SetBetaEquilibrium\", please also set SetBetaEquilibrium::ID_path = </path/to/EOS/slice> (Lorene format) or SetBetaEquilibrium::ID_path = </path/to/EOS/slice>/eos (Compose format).");
    }

    // Basic consistency check 3
    if (not CCTK_EQUALS(ID_path, "")) {
        if (not(CCTK_EQUALS(initial_temperature, "SetBetaEquilibrium")))
            CCTK_PARAMWARN("Since SetBetaEquilibrium::ID_path is (hopefully) set to the path to the beta-equilibrated initial data slice, please also set HydroBase::initial_temperature = \"SetBetaEquilibrium\".");

        if (not(CCTK_EQUALS(initial_Y_e, "SetBetaEquilibrium")))
            CCTK_PARAMWARN("Since SetBetaEquilibrium::ID_path is (hopefully) set to the path to the beta-equilibrated initial data slice, please also set HydroBase::initial_Y_e = \"SetBetaEquilibrium\".");
    }


    /* FIXME: the existence of the following files is anyway checked in
              SetBetaEquilibrium. Should the following be erased?               */
    // Check whether the needed EOS files exist or not
    if (CCTK_EQUALS(ID_format, "Compose")) {
        filepath  = ID_path_string;
        filepath += ".t";
        //if (experimental::filesystem::exists(filepath))  // Not linking on Frontera
        file.open(filepath);
        if (file) {
            file.close();
            CCTK_VINFO("%s found", filepath.c_str());
        } else CCTK_VERROR("%s not found", filepath.c_str());

        filepath  = ID_path_string;
        filepath += ".nb";
        //if (experimental::filesystem::exists(filepath))  // Not linking on Frontera
        file.open(filepath);
        if (file) {
            file.close();
            CCTK_VINFO("%s found", filepath.c_str());
        } else CCTK_VERROR("%s not found", filepath.c_str());


        filepath  = ID_path_string;
        filepath += ".yqb";
        //if (experimental::filesystem::exists(filepath))  // Not linking on Frontera
        file.open(filepath);
        if (file) {
            file.close();
            CCTK_VINFO("%s found", filepath.c_str());
        } else CCTK_VERROR("%s not found", filepath.c_str());


        filepath  = ID_path_string;
        filepath += ".thermo";
        //if (experimental::filesystem::exists(filepath))  // Not linking on Frontera
        file.open(filepath);
        if (file) {
            file.close();
            CCTK_VINFO("%s found", filepath.c_str());
        } else CCTK_VERROR("%s not found", filepath.c_str());
    }


    else if (CCTK_EQUALS(ID_format, "Lorene")) {
        filepath  = ID_path_string;
        filepath += Lorene_slice_filename;
        //if (experimental::filesystem::exists(filepath))  // Not linking on Frontera
        file.open(filepath);
        if (file) {
            file.close();
            CCTK_VINFO("%s found", filepath.c_str());
        } else CCTK_VERROR("%s not found", filepath.c_str());


        filepath  = ID_path_string;
        filepath += ye_filename;
        //if (experimental::filesystem::exists(filepath))  // Not linking on Frontera
        file.open(filepath);
        if (file) {
            file.close();
            CCTK_VINFO("%s found", filepath.c_str());
        } else CCTK_VERROR("%s not found", filepath.c_str());

    }


    else CCTK_ERROR("Either set SetBetaEquilibrium::ID_path = \"Compose\" or SetBetaEquilibrium::ID_path = \"Lorene\"");

    return;
}
