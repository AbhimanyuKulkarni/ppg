#ifndef PPG_DATAGEN_H
#define PPG_DATAGEN_H
#include "setup.h"
// Generate a random sampling matrix
void get_random_sample_flags(size_t N, size_t m, bool *flags);

// Use sampling flags to get compressed samples.
void get_compressed_samples(PPG_FRAC *samples, PPG_Params params,
														bool *phi_flags, PPG_FRAC *compressed);

// Generate simulated cosine wave data from params, into samples.
void generate_sim_data(PPG_Params params, PPG_FRAC *samples);

// Load PhysioNet data into samples.
void get_physionet_data(PPG_FRAC *samples);
void get_physionet_down4_data(PPG_FRAC *samples);
#endif
