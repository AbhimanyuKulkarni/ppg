#ifndef PPG_SETUP_H
#define PPG_SETUP_H

/*
#ifdef PPG_USE_FLOAT64
	#define PPG_FRAC_TYPE double
	#include "float64_ops.h"
#elif defined PPG_USE_FLOAT32
	#define PPG_FRAC_TYPE float
	#include "float32_ops.h"
#elif defined PPG_USE_FXP64
	#define PPG_FRAC_TYPE FxP64
	#include "fxp64_ops.h"
#elif defined PPG_USE_FXP32
	#define PPG_FRAC_TYPE FxP32
	#include "fxp32_ops.h"
#endif
*/

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

	double LASSO_lambda;							// LASSO regularization parameter
	double CD_diff_thresh;						// stopping threshold for coordinate
																		// descent algorithm
	
	char *original_samples_filename;	// where the uncompressed signal is
	char *phi_flags_filename;					// where the compression matrix
																		// flags are
	char *psiT_filename;							// where the DCT matrix is stored
	char *A_mat_filename;							// where the 'A' matrix is stored
	char *A_flip_mat_filename;				// half-flipped 'A'
	char *subsamples_filename;				// where the compressed signal will
																		// be
	char *reconstruction_filename;		// where the reconstructed signal
																		// will go
} PPG_Params;

PPG_Params get_ppg_params(const char *expt_name);
#endif
