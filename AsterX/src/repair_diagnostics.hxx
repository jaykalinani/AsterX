#ifndef ASTERX_REPAIR_DIAGNOSTICS_HXX
#define ASTERX_REPAIR_DIAGNOSTICS_HXX

#include <cctk.h>
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <memory>

namespace AsterX {
enum repair_event {
  cell_atmo, face_atmo, rho_clamp, temp_clamp, ye_clamp, eps_clamp,
  tau_repair, primary_failure, backup_call, backup_failure,
  conservative_recompute, loworder_face, pplim_activation, num_repair_events
};

struct repair_counts {
  unsigned long long value[num_repair_events]{};
};

CCTK_HOST CCTK_DEVICE inline void
count_repair(repair_counts *counts, const repair_event event,
              const unsigned long long amount = 1) {
  if (counts && amount)
    amrex::HostDevice::Atomic::Add(&counts->value[event], amount);
}

// Storage is allocated only for sampled calls. Each invocation owns its
// counters, so concurrent grid calls do not share mutable host state.
class repair_diagnostics {
  std::unique_ptr<amrex::Gpu::DeviceScalar<repair_counts>> storage;

public:
  explicit repair_diagnostics(const bool enabled) {
    if (enabled)
      storage = std::make_unique<amrex::Gpu::DeviceScalar<repair_counts>>(
          repair_counts{});
  }
  repair_counts *data() { return storage ? storage->dataPtr() : nullptr; }

  void report(const cGH *cctkGH, const char *where, const int dir = -1) const {
    if (!storage)
      return;
    amrex::Gpu::streamSynchronize();
    const auto counts = storage->dataValue();
    CCTK_VINFO("EOS repairs: %s rank=%d iteration=%d dir=%d "
               "cell_atmo=%llu face_atmo=%llu rho_clamp=%llu T_clamp=%llu "
               "Ye_clamp=%llu eps_clamp=%llu tau_repair=%llu "
               "primary_fail=%llu backup_call=%llu backup_fail=%llu "
               "cons_recompute=%llu loworder_face=%llu pplim=%llu",
               where, CCTK_MyProc(cctkGH), int(cctkGH->cctk_iteration), dir,
               counts.value[cell_atmo], counts.value[face_atmo],
               counts.value[rho_clamp], counts.value[temp_clamp],
               counts.value[ye_clamp], counts.value[eps_clamp],
               counts.value[tau_repair], counts.value[primary_failure],
               counts.value[backup_call], counts.value[backup_failure],
               counts.value[conservative_recompute], counts.value[loworder_face],
               counts.value[pplim_activation]);
  }
};
} // namespace AsterX
#endif
