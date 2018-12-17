#ifndef PPG_ALGORITHM_H
#define PPG_ALGORITHM_H

#include "setup.h"
#include <stdlib.h>
#include <stdbool.h>

// Main algorithm which should be used to reconstruct an entire signal.
void ppg(double *Y, bool *phi_flags, PPG_Params params,
				double *tr, double *Xr);

// Parts of the overall PPG reconstruction algorithm, exposed for
// testing.

// LASSO solver based on coordinate descent.
void cd_lasso(double *y, double *A, size_t M, size_t N,
							double lambda, double tol, double *x_hat);

// Get DCT matrix
void setup_psiT(double *psiT, size_t N);

// Form sampled DCT matrix (A matrix)
void get_selected_rows(double *mat, bool *flags,
											size_t N, size_t m, double *out);

// Half-flip to handle overlapping windows which don't come from the
// same model
void half_flip(bool *flags, size_t N, size_t M, bool *flipped);

// Estimate dominant frequency of ppg signal, uses RMS crossing rate
// subject to max frequency constraint
// Result will be in Hz, multiply by 60 to get bpm
// This is for consistency with the max_f passed in, which should also
// be in Hz.
double get_freq(double *X, double *t,
								size_t N, double max_f);

// Get Pearson correlation with ground-truth signal
double corrcoef(double *a, double *b, size_t N);
#endif
