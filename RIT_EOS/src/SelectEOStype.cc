#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>
#include <RIT_EOS.hh>



// Definition of namespace "RIT_EOS_function_callers" (declared in RIT_EOS.hh)
namespace RIT_EOS_function_callers {
    void (*press_eps_from_rho_temp_ye_caller)(
        const CCTK_REAL &,
        const CCTK_REAL &,
        const CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &);

    void (*press_temp_from_rho_ye_eps_caller)(
        const CCTK_REAL &,
        const CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &);

    void (*dpdrhoe_temp_from_rho_ye_eps_caller)(
        const CCTK_REAL &, 
        const CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &);

    void (*dpderho_temp_from_rho_ye_eps_caller)(
        const CCTK_REAL &, 
        const CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &);

    void (*dpdrhoe_dpderho_temp_from_rho_ye_eps_caller)(
        const CCTK_REAL &, 
        const CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &,
              CCTK_REAL &);
}
    
using namespace RIT_EOS_function_callers;





extern "C" void SelectEOStype(CCTK_ARGUMENTS) {
    DECLARE_CCTK_PARAMETERS;

    switch (EOS_type) {
        // Polytropic EOS (this includes the ideal-fluid EOS as a special case)
        case 1:
            press_eps_from_rho_temp_ye_caller           = press_eps_from_rho_temp_ye_polyEOS;
            press_temp_from_rho_ye_eps_caller           = press_temp_from_rho_ye_eps_polyEOS;
            dpdrhoe_temp_from_rho_ye_eps_caller         = dpdrhoe_temp_from_rho_ye_eps_polyEOS;
            dpderho_temp_from_rho_ye_eps_caller         = dpderho_temp_from_rho_ye_eps_polyEOS;
            dpdrhoe_dpderho_temp_from_rho_ye_eps_caller = dpdrhoe_dpderho_temp_from_rho_ye_eps_polyEOS;
            CCTK_INFO("Polytropic EOS selected");
            break;

        // Finite-temperature, microphysical, tabulated EOS
        case 2:
            press_eps_from_rho_temp_ye_caller           = press_eps_from_rho_temp_ye_nucEOS;
            press_temp_from_rho_ye_eps_caller           = press_temp_from_rho_ye_eps_nucEOS;
            dpdrhoe_temp_from_rho_ye_eps_caller         = dpdrhoe_temp_from_rho_ye_eps_nucEOS;
            dpderho_temp_from_rho_ye_eps_caller         = dpderho_temp_from_rho_ye_eps_nucEOS;
            dpdrhoe_dpderho_temp_from_rho_ye_eps_caller = dpdrhoe_dpderho_temp_from_rho_ye_eps_nucEOS;
            CCTK_INFO("Finite-temperature, microphysical, tabulated EOS selected");
            break;

        default:
            CCTK_VERROR("Parameter EOS_type = %d is not supported.", EOS_type);
    }

    return;
}
