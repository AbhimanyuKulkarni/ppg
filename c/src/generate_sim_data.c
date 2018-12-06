/*
 * generate_sim_data.c
 *
 * Generate simulated data for use by PPG. Useful for when files are not
 * available.
 *
 * The simulated data will be an amplitude-modulated and noisy cosine
 * wave of base frequency 1 Hz. Gaussian noise is too much work to
 * write, so we'll go with uniform(0, 2) noise instead.
 */

#include "setup.h"
#include <time.h>
#include <stdlib.h>
#include <math.h>

void generate_sim_data(PPG_Params params, double *samples) {
	srand(time(NULL));

	size_t smoothing_factor = 0.9;				// make AM low-pass
	double amp_A = 0.1;										// for AM
	double sigma0 = 0.15;									// adding random noise
	double b0 = (((double) rand()) / RAND_MAX) * 2 * M_PI;
																				// random initial phase
	
	double A = 0;
	double u0 = 0;
	double t = 0;
	for (size_t i = 0; i < params.N0; ++i) {
		A *= smoothing_factor;
		A += (1 - smoothing_factor) * 
					(1.0 + (((double) rand()) / RAND_MAX - 0.5) * amp_A);
		u0 = 2 * (((double) rand()) / RAND_MAX);
		samples[i] = A * (cos(2 * M_PI * t - b0) + u0);
		t += params.T0;
	}
}
