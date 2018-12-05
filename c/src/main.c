/*
 * Runs the reconstruction algorithm from ppg_algorithm.c on some file.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

#include "setup.h"
#include "ppg_algorithm.h"

double max(double a, double b) {
	return (a > b ? a : b) ;
}


void get_random_sample_flags_from_file(char *filepath, size_t N, bool *flags) {
	FILE *fp = fopen(filepath, "r");
	if (fp) {
		int val;
		for (int i = 0; i < N; ++i) {
			// I'm afraid of int value clobbering many bools in one shot
			// somehow.
			// I don't want to bother checking that though, so val it is.
			fscanf(fp, "%d\n", &val);
			flags[i] = (bool) val;
		}
	}
	fclose(fp);
}

double corrcoef(double *a, double *b, size_t N) {
	double mean_a, mean_b;
	mean_a = mean_b = 0.0;
	for (size_t i = 0; i < N; ++i) {
		mean_a += (1.0 / N) * a[i];
		mean_b += (1.0 / N) * b[i];
	}

	double norm_a = 0.0, norm_b = 0.0;
	for (size_t i = 0; i < N; ++i) {
		norm_a += pow(a[i] - mean_a, 2);
		norm_b += pow(b[i] - mean_b, 2);
	}
	norm_a = sqrt(norm_a);
	norm_b = sqrt(norm_b);

	// Don't let (1/std) blow up; better to just let it die.
	if (norm_a < 1.0e-4) norm_a = 1.0e10;
	if (norm_b < 1.0e-4) norm_b = 1.0e10;

	double corr = 0.0;
	for (size_t i = 0; i < N; ++i) {
		corr += ((a[i] - mean_a) * (b[i] - mean_b));
	}
	corr /= (norm_a * norm_b);
	return corr;
}


/*
 * argv[0]
 * argv[1]: <expt_name> (physionet or sim)
 */
int main(int argc, char **argv) {
	if (argc < 2) {
		fprintf(stderr, "Not enough arguments!\n");
		fprintf(stderr, "argv[1]: <expt_name> (physionet or sim)\n");
	}
	
	srand(0);

	// from experiment_setup import *
	PPG_Params params = get_ppg_params(argv[1]);
	printf("f0: %ld\nT0: %.3lf\n", params.f0, params.T0);
	printf("N_window: %ld\nM: %ld\n", params.N_window, params.M);
	
	// original_data = np.loadtxt(original_samples_filename,
	//														delimiter=',')
	double *t0, *ts, *X0, *Y, *Xr;
	t0 = malloc(sizeof(double) * params.N0);
	X0 = malloc(sizeof(double) * params.N0);
	ts = malloc(sizeof(double) * params.N0 / params.CF);
	Y  = malloc(sizeof(double) * params.N0 / params.CF);
	Xr = malloc(sizeof(double) * params.N0);

	size_t n = 0;
	FILE* fp_samples = fopen(params.original_samples_filename, "r");
	if (!fp_samples) {
		fprintf(stderr, "Couldn't open file %s!\n",
						params.original_samples_filename);
		exit(EXIT_FAILURE);
	} else {
		while (fscanf(fp_samples, "%lf,%lf\n", &t0[n], &X0[n]) == 2) {
			n++;
			if (n == params.N0) break;
		}
		fclose(fp_samples);
	}
	printf("Read %ld samples from %s.\n",
				n, params.original_samples_filename);

	// subsampled_data = np.loadtxt(params.subsamples_filename,
	//															delimiter=',')
	n = 0;
	fp_samples = fopen(params.subsamples_filename, "r");
	if (!fp_samples) {
		fprintf(stderr, "Couldn't open file %s!\n",
						params.subsamples_filename);
	} else {
		while (fscanf(fp_samples, "%lf,%lf\n", &ts[n], &Y[n]) == 2) {
			n++;
			if (n == params.N0 / params.CF) break;
		}
		fclose(fp_samples);
	}
	printf("Read %ld samples from %s.\n", n, params.subsamples_filename);

	// phi_flags = np.loadtxt(phi_flags_filename, delimiter=',')
	bool *phi_flags = malloc(sizeof(bool) * params.N_window);
	get_random_sample_flags_from_file(params.phi_flags_filename,
																		params.N_window, phi_flags);

	// TODO: period estimation
	
	// call ppg algorithm
	ppg(Y, phi_flags, params, Xr);
	
	printf("Correlation coefficient: %.4lf\n",
					corrcoef(X0, Xr, params.N0));

	// np.savetxt(reconstruction_filename, out_data, delimiter=',')
	printf("Writing to file %s...\n", params.reconstruction_filename);
	FILE *fp_write = fopen(params.reconstruction_filename, "w");
	if (!fp_write) {
		fprintf(stderr, "Couldn't open file %s for writing!\n", argv[2]);
	} else {
		for (size_t i = 0; i < params.N0; ++i) {
			fprintf(fp_write, "%f,%f\n", t0[i], Xr[i]);
		}
	}
	fclose(fp_write);

	free(t0);
	free(ts);
	free(X0);
	free(Y);
	free(Xr);
	free(phi_flags);
	return 0;
}
