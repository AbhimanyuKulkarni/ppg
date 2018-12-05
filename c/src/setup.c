#include "setup.h"
#include <string.h>

size_t nearest_multiple(size_t N, size_t d) {
	return N - (N % d);
}

PPG_Params get_ppg_params(const char *expt_name) {
	PPG_Params params;
	params.MAX_SAMPLES = 10000;
	if (strcmp(expt_name, "physionet") == 0) {
		params.T = 60;
		params.f0 = 256;
		params.T0 = 1.0 / params.f0;
		params.T_window = 60;
		params.N_window = nearest_multiple(params.T_window * params.f0, 4);
		params.CF = 32;
		params.M = params.N_window / params.CF;
		params.N0 = nearest_multiple(params.T * params.f0,
																 params.N_window / 2);
		
		params.original_samples_filename = "data/physionet/samples.csv";
		params.phi_flags_filename
						= "data/physionet/phi_flags_15360_480.csv";
		params.subsamples_filename
						= "data/physionet/subsamples_15360_480.csv";
		params.reconstruction_filename
						= "data/physionet/reconstructed_c.csv";
	} else if (strcmp(expt_name, "sim") == 0) {
		params.T = 600;
		params.f0 = 4;
		params.T0 = 1.0 / params.f0;
		params.T_window = 60;
		params.N_window = nearest_multiple(params.T_window * params.f0, 4);
		params.CF = 12;
		params.M = params.N_window / params.CF;
		params.N0 = nearest_multiple(params.T * params.f0,
																 params.N_window / 2);
		
		params.original_samples_filename = "data/sim/samples.csv";
		params.phi_flags_filename
						= "data/sim/phi_flags_240_20.csv";
		params.subsamples_filename
						= "data/sim/subsamples_240_20.csv";
		params.reconstruction_filename
						= "data/sim/reconstructed_c.csv";
	}

	params.LASSO_lambda = 1.0e-2;
	params.CD_diff_thresh = 1.0e-2;
	return params;
}
