#ifndef PPG_DATAGEN_H
#define PPG_DATAGEN_H
#include "setup.h"
// Generate a random sampling matrix
void get_random_sample_flags(size_t N, size_t m, bool *flags);

// Use sampling flags to get compressed samples.
void get_compressed_samples(double *samples, PPG_Params params,
														bool *phi_flags, double *compressed);

// Generate simulated cosine wave data from params, into samples.
void generate_sim_data(PPG_Params params, double *samples);

// Load PhysioNet data into samples.
void get_physionet_data(double *samples);
void get_physionet_down4_data(double *samples);
#endif
