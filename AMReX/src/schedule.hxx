#ifndef SCHEDULE_HXX
#define SCHEDULE_HXX

#include <cctk.h>
#include <cctk_Schedule.h>

namespace AMReX {

void setup_cctkGH(cGH *restrict cctkGH);
void enter_global_mode(cGH *restrict cctkGH);
void leave_global_mode(cGH *restrict cctkGH);
void enter_level_mode(cGH *restrict cctkGH);
void leave_level_mode(cGH *restrict cctkGH);
void enter_local_mode(cGH *restrict cctkGH, const MFIter &mfi);
void leave_local_mode(cGH *restrict cctkGH);

int CallFunction(void *function, cFunctionData *attribute, void *data);
} // namespace AMReX

#endif // #ifndef SCHEDULE_HXX
