#include "timer.hxx"

#include <cassert>

namespace CarpetX {

Timer::Timer(const string &name)
    : name(name), handle(CCTK_TimerCreate(name.c_str())) {}

void Timer::print() const {
  const int is_running = CCTK_TimerIsRunningI(handle);
  assert(is_running >= 0);
  // We need to stop timers before printing them; running timers cannot be
  // printed
  if (is_running)
    CCTK_TimerStopI(handle);
  const int num_clocks = CCTK_NumClocks();
  for (int clock = 0; clock < num_clocks; ++clock)
    CCTK_TimerPrintDataI(handle, clock);
  if (is_running)
    CCTK_TimerStartI(handle);
}

Interval::Interval(const Timer &timer) : timer(timer) {
  CCTK_TimerStartI(timer.handle);
#ifdef __CUDACC__
  range = nvtxRangeStartA(timer.name.c_str());
#endif
}

Interval::~Interval() {
  CCTK_TimerStopI(timer.handle);
#ifdef __CUDACC__
  nvtxRangeEnd(range);
#endif
}

} // namespace CarpetX
