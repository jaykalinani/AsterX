#include <cctk.h>
#include <cctk_Arguments.h>

#include <Kokkos_Core.hpp>

extern "C" int Kokkos_Startup() {
  int argc = 0;
  char *argv[] = {nullptr};
  Kokkos::initialize(argc, argv);
  return 0;
}

extern "C" int Kokkos_Shutdown() {
  Kokkos::finalize();
  return 0;
}
