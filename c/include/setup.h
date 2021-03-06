#ifndef PPG_SETUP_H
#define PPG_SETUP_H

#if defined(PPG_USE_FXP64) || defined(PPG_USE_FXP32) \
		|| defined(PPG_USE_FXP16) || defined(PPG_USE_FXP8)
	// use fixed-point
	#define PPG_FRAC FxP
	#include "fxp_ops.h"
#elif defined(PPG_USE_FLOAT32)
	#define PPG_FRAC float
	#include "float_ops.h"
#elif defined(PPG_USE_FLOAT64)
	#define PPG_FRAC double
	#include "float_ops.h"
#else	// default
	#define PPG_FRAC double
	#include "float_ops.h"
#endif

#include <stddef.h>

typedef struct PPG_Params {
	size_t MAX_SAMPLES;	// maximum number of samples we can handle

	size_t T;						// total length of signal, in seconds
	size_t N0;					// #samples in whole signal

	size_t f0;					// sample rate of original signal, in Hz
	PPG_FRAC T0;				// 1 / f0, so in seconds

	size_t T_window;		// Window length for reconstruction algorithm,
											// in seconds
	size_t N_window;		// #samples in that window

	unsigned int CF;		// compression factor
	size_t M;						// compressed window length = N_window / CF

	PPG_FRAC LASSO_lambda;						// LASSO regularization parameter
	PPG_FRAC CD_diff_thresh;					// stopping threshold for coordinate
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
