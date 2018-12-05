#ifndef PPG_ALGORITHM_H
#define PPG_ALGORITHM_H

#include <stdlib.h>
#include <stdbool.h>
#include "setup.h"

// Main algorithm which should be used to reconstruct an entire signal.
void ppg(double *Y, bool *phi_flags, PPG_Params params, double *Xr);
#endif
