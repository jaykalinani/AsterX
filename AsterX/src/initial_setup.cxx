#include <initial_setup.hxx>

namespace AsterX {

using namespace EOSX;

enum class eos_id {Polytropic, PWPolytropic}
enum class eos_evol { IdealGas, Hybrid, Tabulated };

// initial data EOS
Vector<eos_polytrope, 1> eos_poly;

// evolution EOS
Vector<eos_idealgas, 1> eos_ig;
Vector<eos_tabulated3d, 1> eos_tab3d;

extern "C" void AsterX_Initialize_EOSID(CCTK_ARGUMENTS) {

  eos_id eos_id_type;

  if (CCTK_EQUALS(id_eos_name, "Polytropic")) {
    eos_id_type = eos_id::Polytropic;
  } else if (CCTK_EQUALS(id_eos_name, "PWPolytropic")) {
    eos_id_type = eos_id::PWPolytropic;
  } else {
    CCTK_ERROR("Unknown value for parameter \"initial_data_eos\"");
  }

  switch (eos_id_type) {
    case eos_id::Polytropic: {
      eos_poly.init(poly_gamma, poly_k, rho_max);
    }
    case eos_id::PWPolytropic: {
      CCTK_ERROR("Piecewise Polytrope EOS is not supported yet!");
    }
    default:
      assert(0);
  }
}

extern "C" void AsterX_Initialize_EOS(CCTK_ARGUMENTS) {

  eos_evol eos_evol_type;
  if (CCTK_EQUALS(evol_eos_name, "IdealGas")) {
    eos_evol_type = eos_evol::IdealGas;
  } else if (CCTK_EQUALS(evol_eos_name, "Hybrid")) {
    eos_evol_type = eos_evol::Hybrid;
  } else if (CCTK_EQUALS(evol_eos_name, "Tabulated3d")) {
    eos_evol_type = eos_evol::Tabulated;
  } else {
    CCTK_ERROR("Unknown value for parameter \"evolution_eos\"");
  }

  switch (eos_evol_type) {
  case eos_evol::IdealGas: {
    eos_ig.init(gl_gamma, particle_mass, rgeps, rgrho, rgye);
    break;
  }
  case eos_evol::Hybrid: {
    CCTK_ERROR("Hybrid EOS is not yet supported");
    break;
  }
  case eos_evol::Tabulated: {
    eos_tab3d.init(rgeps, rgrho, rgye);
    const string eos_filename = EOSTable_filename;
    eos_tab3d.read_eos_table(eos_filename);
    break;
  }
  default:
    assert(0);
  }

}

} // namespace AsterX

#endif