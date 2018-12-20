#ifndef PPG_ALGORITHM_H
#define PPG_ALGORITHM_H

#include "setup.h"
#include <stdlib.h>
#include <stdbool.h>

// Main algorithm which should be used to reconstruct an entire signal.
void ppg(PPG_FRAC *Y, bool *phi_flags, PPG_Params params,
				PPG_FRAC *tr, PPG_FRAC *Xr);

// Parts of the overall PPG reconstruction algorithm, exposed for
// testing.

// Standardization x <- (x - min) / (max - min), with further constraint
// on range from val_range. Final range [-1/M, 1/M] where 
// M = max(1, N/val_range)
void standardize(PPG_FRAC *data, size_t N, size_t val_range);

// LASSO solver based on coordinate descent.
void cd_lasso(PPG_FRAC *y, PPG_FRAC *A, size_t M, size_t N,
							PPG_FRAC lambda, PPG_FRAC tol, PPG_FRAC *x_hat);

// Get DCT matrix
void setup_psiT(PPG_FRAC *psiT, size_t N);

// Form sampled DCT matrix (A matrix)
void get_selected_rows(PPG_FRAC *mat, bool *flags,
											size_t N, size_t m, PPG_FRAC *out);

// Half-flip to handle overlapping windows which don't come from the
// same model
void half_flip(bool *flags, size_t N, size_t M, bool *flipped);

// Estimate dominant frequency of ppg signal, uses RMS crossing rate
// subject to max frequency constraint
// Result will be in Hz, multiply by 60 to get bpm
// This is for consistency with the max_f passed in, which should also
// be in Hz.
double get_freq(PPG_FRAC *X, PPG_FRAC *t,
								size_t N, double max_f);

// Get Pearson correlation with ground-truth signal
double corrcoef(PPG_FRAC *a, PPG_FRAC *b, size_t N);
#endif
