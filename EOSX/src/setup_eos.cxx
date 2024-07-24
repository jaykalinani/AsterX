#include <cctk.h>
#include <cctk_Arguments.h>
#include <cctk_Parameters.h>

#include <AMReX.H>

#include <setup_eos.hxx>


namespace EOSX {

using namespace amrex;

enum class eos_1param {Polytropic, PWPolytropic};
enum class eos_3param { IdealGas, Hybrid, Tabulated };

// initial data EOS
// AMREX_GPU_MANAGED eos_1p_polytrope *eos_1p_poly = nullptr;

// evolution EOS
// AMREX_GPU_MANAGED eos_3p_idealgas    *eos_3p_ig    = nullptr;
// AMREX_GPU_MANAGED eos_3p_hybrid      *eos_3p_hyb    = nullptr;
// AMREX_GPU_MANAGED eos_3p_tabulated3d *eos_3p_tab3d = nullptr;


extern "C" void EOSX_Setup_EOSID(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  eos_1param eos_1p_type;
  if (CCTK_EQUALS(initial_data_eos, "Polytropic")) {
    eos_1p_type = eos_1param::Polytropic;
  } else if (CCTK_EQUALS(initial_data_eos, "PWPolytropic")) {
    eos_1p_type = eos_1param::PWPolytropic;
  } else {
    CCTK_ERROR("Unknown value for parameter \"initial_data_eos\"");
  }

  switch (eos_1p_type) {
    case eos_1param::Polytropic: {
      CCTK_INFO("Setting initial data EOS to Polytropic");
      eos_1p_poly = (eos_1p_polytrope*)The_Managed_Arena()->alloc(sizeof *eos_1p_poly);
      new (eos_1p_poly) eos_1p_polytrope;
      assert(eos_1p_poly);
      eos_1p_poly->init(poly_gamma, poly_k, rho_max);
      break;
    }
    case eos_1param::PWPolytropic: {
      CCTK_ERROR("Piecewise Polytrope EOS is not supported yet!");
      break;
    }
    default:
      assert(0);
  }
}

extern "C" void EOSX_Setup_EOS(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  eos_3param eos_3p_type;
  eos_3p::range rgeps(eps_min, eps_max), rgrho(rho_min, rho_max),
      rgye(ye_min, ye_max);

  if (CCTK_EQUALS(evolution_eos, "IdealGas")) {
    eos_3p_type = eos_3param::IdealGas;
  } else if (CCTK_EQUALS(evolution_eos, "Hybrid")) {
    eos_3p_type = eos_3param::Hybrid;
  } else if (CCTK_EQUALS(evolution_eos, "Tabulated3d")) {
    eos_3p_type = eos_3param::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_3p_type) {
  case eos_3param::IdealGas: {
    CCTK_INFO("Setting evolution EOS to Ideal Gas");
    eos_3p_ig = (eos_3p_idealgas*)The_Managed_Arena()->alloc(sizeof *eos_3p_ig);
    new (eos_3p_ig) eos_3p_idealgas;
    assert(eos_3p_ig);
    eos_3p_ig->init(gl_gamma, particle_mass, rgeps, rgrho, rgye);
    break;
  }
  case eos_3param::Hybrid: {
    CCTK_INFO("Setting evolution EOS to Hybrid");
    eos_3p_hyb = (eos_3p_hybrid*)The_Managed_Arena()->alloc(sizeof *eos_3p_hyb);
    new (eos_3p_hyb) eos_3p_hybrid(eos_1p_poly, gamma_th, rgeps, rgrho, rgye);
    assert(eos_3p_hyb);
    break;
  }
  case eos_3param::Tabulated: {
    CCTK_INFO("Setting evolution EOS to Tabulated3D");
    const string eos_filename = EOSTable_filename;
    eos_3p_tab3d = (eos_3p_tabulated3d*)The_Managed_Arena()->alloc(sizeof *eos_3p_tab3d);
    new (eos_3p_tab3d) eos_3p_tabulated3d;
    assert(eos_3p_tab3d);
    eos_3p_tab3d->init(rgeps, rgrho, rgye);
    eos_3p_tab3d->read_eos_table(eos_filename);
    break;
  }
  default:
    assert(0);
  }

}

} // namespace EOSX

