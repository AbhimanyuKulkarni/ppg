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
#include "datagen.h"

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

/*
 * argv[0]
 * argv[1]: <expt_name> (physionet, physionet_down4, or sim)
 * argv[2]: <mode> (file or gen)
 */
int main(int argc, char **argv) {
	if (argc < 3) {
		fprintf(stderr, "Not enough arguments!\n");
		fprintf(stderr, "argv[1]: <expt_name> (physionet, physionet_down4, or sim)\n");
		fprintf(stderr, "argv[2]: <mode> (file or gen)\n");
		exit(1);
	}
	
	srand(time(NULL));

	// from experiment_setup import *
	PPG_Params params = get_ppg_params(argv[1]);
	bool use_datagen = true;
	if (strcmp(argv[2], "file") == 0) {
		use_datagen = false;
	} else if (strcmp(argv[2], "gen") == 0) {
		use_datagen = true;
	}

	printf("f0: %ld\nT0: %.3lf\n", params.f0, FRAC_TO_FLOAT64(params.T0));
	printf("N_window: %ld\nM: %ld\n", params.N_window, params.M);
	
	// original_data = np.loadtxt(original_samples_filename,
	//														delimiter=',')
	PPG_FRAC *t0, *ts, *X0, *Y, *tr, *Xr;
	t0 = malloc(sizeof(PPG_FRAC) * params.N0);
	X0 = malloc(sizeof(PPG_FRAC) * params.N0);
	ts = malloc(sizeof(PPG_FRAC) * params.N0 / params.CF);
	Y  = malloc(sizeof(PPG_FRAC) * params.N0 / params.CF);
	tr = malloc(sizeof(PPG_FRAC) * params.N0);
	Xr = malloc(sizeof(PPG_FRAC) * params.N0);
	bool *phi_flags = malloc(sizeof(bool) * params.N_window);

	if (!use_datagen) {
		size_t n = 0;
		FILE* fp_samples = fopen(params.original_samples_filename, "r");
		if (!fp_samples) {
			fprintf(stderr, "Couldn't open file %s!\n",
							params.original_samples_filename);
			exit(EXIT_FAILURE);
		} else {
			double t, x;
			while (fscanf(fp_samples, "%lf,%lf\n", &t, &x) == 2) {
				t0[n] = FLOAT64_TO_FRAC(t);
				X0[n] = FLOAT64_TO_FRAC(x);
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
			double t, y;
			while (fscanf(fp_samples, "%lf,%lf\n", &t, &y) == 2) {
				ts[n] = FLOAT64_TO_FRAC(t);
				Y[n] = FLOAT64_TO_FRAC(y);
				n++;
				if (n == params.N0 / params.CF) break;
			}
			fclose(fp_samples);
		}
		printf("Read %ld samples from %s.\n", n,
					params.subsamples_filename);

		// phi_flags = np.loadtxt(phi_flags_filename, delimiter=',')
		get_random_sample_flags_from_file(params.phi_flags_filename,
																			params.N_window, phi_flags);

	} else {
		if (strcmp(argv[1], "sim") == 0) {
			generate_sim_data(params, X0);
		} else if (strcmp(argv[1], "physionet") == 0) {
			get_physionet_data(X0);
		} else if (strcmp(argv[1], "physionet_down4") == 0) {
			get_physionet_down4_data(X0);
		}

		get_random_sample_flags(params.N_window, params.M, phi_flags);
		get_compressed_samples(X0, params, phi_flags, Y);
	}

	// call ppg algorithm
	ppg(Y, phi_flags, params, tr, Xr);

	printf("Correlation coefficient: %.4lf\n",
					corrcoef(X0, Xr, params.N0));

	// np.savetxt(reconstruction_filename, out_data, delimiter=',')
	printf("Writing to file %s...\n", params.reconstruction_filename);
	FILE *fp_write = fopen(params.reconstruction_filename, "w");
	if (!fp_write) {
		fprintf(stderr, "Couldn't open file %s for writing!\n", argv[2]);
	} else {
		for (size_t i = 0; i < params.N0; ++i) {
			fprintf(fp_write, "%f,%f\n",
							FRAC_TO_FLOAT64(tr[i]), FRAC_TO_FLOAT64(Xr[i]));
		}
	}
	fclose(fp_write);

	free(t0);
	free(ts);
	free(X0);
	free(Y);
	free(tr);
	free(Xr);
	free(phi_flags);
	return 0;
}
