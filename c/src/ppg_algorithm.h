#ifndef PPG_ALGORITHM_H
#define PPG_ALGORITHM_H

#include <stdlib.h>
#include <stdbool.h>
#include "setup.h"

// Main algorithm which should be used to reconstruct an entire signal.
void ppg(double *Y, bool *phi_flags, PPG_Params params,
				double *tr, double *Xr);

// Estimate dominant frequency of ppg signal, uses RMS crossing rate
// subject to max frequency constraint
// Result will be in Hz, multiply by 60 to get bpm
// This is for consistency with the max_f passed in, which should also
// be in Hz.
double get_freq(double *X, double *t, size_t N, double max_f);
#endif
