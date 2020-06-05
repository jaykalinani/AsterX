#ifndef TIMER_HXX
#define TIMER_HXX

#include <cctk.h>

#include <string>

namespace CarpetX {
using namespace std;

class Interval;

class Timer {
  friend class Interval;

  string name;
  int handle;

public:
  Timer() = delete;
  Timer(const string &name);
};

class Interval {
  const Timer &timer;

public:
  Interval() = delete;
  Interval(const Timer &timer);
  ~Interval();
};

} // namespace CarpetX

#endif // #ifndef TIMER_HXX
