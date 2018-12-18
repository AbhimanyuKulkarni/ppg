#include "setup.h"
#include <string.h>

size_t nearest_multiple(size_t N, size_t d) {
	return N - (N % d);
}

PPG_Params get_ppg_params(const char *expt_name) {
	PPG_Params params;
	params.MAX_SAMPLES = 16384;
	if (strcmp(expt_name, "physionet") == 0) {
		params.T = 60;
		params.f0 = 256;
		params.T_window = 60;
		params.CF = 32;
		
		params.original_samples_filename = "data/physionet/samples.csv";
		params.phi_flags_filename
						= "data/physionet/phi_flags_15360_480.csv";
		params.subsamples_filename
						= "data/physionet/subsamples_15360_480.csv";
		params.reconstruction_filename
						= "data/physionet/reconstructed_c.csv";
		params.psiT_filename
						= "data/physionet/psiT_15360.csv";
		params.A_mat_filename
						= "data/physionet/A_15360_480.csv";
		params.A_flip_mat_filename
						= "data/physionet/A_flip_15360_480.csv";
	} else if (strcmp(expt_name, "physionet_down4") == 0) {
		params.T = 60;
		params.f0 = 64;
		params.T_window = 60;
		params.CF = 32;
		
		params.original_samples_filename = "data/physionet/samples_down4.csv";
		params.phi_flags_filename
						= "data/physionet/phi_flags_3840_120.csv";
		params.subsamples_filename
						= "data/physionet/subsamples_3840_120.csv";
		params.reconstruction_filename
						= "data/physionet/down4_reconstructed_c.csv";
		params.psiT_filename
						= "data/physionet/psiT_3840.csv";
		params.A_mat_filename
						= "data/physionet/A_3840_120.csv";
		params.A_flip_mat_filename
						= "data/physionet/A_flip_3840_120.csv";
	} else if (strcmp(expt_name, "sim") == 0) {
		params.T = 600;
		params.f0 = 4;
		params.T_window = 60;
		params.CF = 12;
		
		params.original_samples_filename = "data/sim/samples.csv";
		params.phi_flags_filename
						= "data/sim/phi_flags_240_20.csv";
		params.subsamples_filename
						= "data/sim/subsamples_240_20.csv";
		params.reconstruction_filename
						= "data/sim/reconstructed_c.csv";
		params.psiT_filename
						= "data/sim/psiT_240.csv";
		params.A_mat_filename
						= "data/sim/A_240_20.csv";
		params.A_flip_mat_filename
						= "data/sim/A_flip_240_20.csv";
	}

	params.T0 = FLOAT64_TO_FRAC(1.0 / params.f0);
	params.N_window = nearest_multiple(params.T_window * params.f0, 4);
	params.M = params.N_window / params.CF;
	params.N0 = nearest_multiple(params.T * params.f0,
															 params.N_window / 2);
	params.LASSO_lambda = FLOAT64_TO_FRAC(1.0e-2);
	params.CD_diff_thresh = FLOAT64_TO_FRAC(1.0e-2);
	return params;
}
