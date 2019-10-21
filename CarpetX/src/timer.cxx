#include "timer.hxx"

namespace CarpetX {

Timer::Timer(const string &name)
    : name(name), handle(CCTK_TimerCreate(name.c_str())) {}

Interval::Interval(const Timer &timer) : timer(timer) {
  CCTK_TimerStartI(timer.handle);
}

Interval::~Interval() { CCTK_TimerStopI(timer.handle); }

} // namespace CarpetX
