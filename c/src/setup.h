#ifndef PPG_SETUP_H
#define PPG_SETUP_H

#include <stddef.h>

typedef struct PPG_Params {
	size_t MAX_SAMPLES;	// maximum number of samples we can handle

	size_t T;						// total length of signal, in seconds
	size_t N0;					// #samples in whole signal

	size_t f0;					// sample rate of original signal, in Hz
	double T0;					// 1 / f0, so in seconds

	size_t T_window;		// Window length for reconstruction algorithm,
											// in seconds
	size_t N_window;		// #samples in that window

	unsigned int CF;		// compression factor
	size_t M;						// compressed window length = N_window / CF

	double LASSO_lambda;// regularization parameter for LASSO
	double CD_diff_thresh;	// stopping threshold for coordinate descent
													// algorithm
	
	char *original_samples_filename;	// where the uncompressed signal is
	char *phi_flags_filename;					// where the compression matrix
																		// flags are
	char *subsamples_filename;				// where the compressed signal will
																		// be
	char *reconstruction_filename;		// where the reconstructed signal
																		// will go
} PPG_Params;

PPG_Params get_ppg_params(const char *expt_name);
#endif
