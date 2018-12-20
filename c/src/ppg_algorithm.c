#include "ppg_algorithm.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

////////////////////////////////////////////////////////////////////////
// Preprocessing
void standardize(PPG_FRAC *data, size_t N, size_t val_range) {
	PPG_FRAC min, max, data_range;
	min = data[0];
	max = data[0];
	const PPG_FRAC ONE = INT64_TO_FRAC(1);
	const PPG_FRAC TWO = INT64_TO_FRAC(2);
	for (size_t i = 0; i < N; ++i) {
		min = MIN(min, data[i]);
		max = MAX(max, data[i]);
	}

	// scale to [-1 / M, 1 / M], where M = max(1, N / val_range)
	data_range = SUB(max, min);
	PPG_FRAC final_range = MAX(ONE, INT64_TO_FRAC(N/val_range));
	for (size_t i = 0; i < N; ++i) {
		data[i] = SUB(MULT(TWO, DIV(SUB(data[i], min), data_range)), ONE);
		data[i] = DIV(data[i], final_range);
	}
}

// Matrix math
static void dot(PPG_FRAC *A, PPG_FRAC *B, size_t M, size_t N, size_t P,
								PPG_FRAC *C) {
	// C = AB, everything is row-major
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < P; ++j) {
			// C[i * P + j] = 0;
			C[i * P + j] = INT64_TO_FRAC(0);
			for (size_t k = 0; k < N; ++k) {
				// C[i * P + j] += A[i * N + k] * B[k * P + j];
				C[i * P + j] = ADD(C[i * P + j],
													 MULT(A[i * N + k], B[k * P + j]));
			}
		}
	}
}

// LASSO solvers
void cd_lasso(PPG_FRAC *y, PPG_FRAC *A, size_t M, size_t N,
							PPG_FRAC lambda, PPG_FRAC tol, PPG_FRAC *x_hat) {
	PPG_FRAC *A_norm2 = malloc(sizeof(PPG_FRAC) * N);
	PPG_FRAC *r = malloc(sizeof(PPG_FRAC) * N); // residual
	const PPG_FRAC ZERO = INT64_TO_FRAC(0);
	const PPG_FRAC ONE = INT64_TO_FRAC(1);
	const PPG_FRAC MINUS_ONE = INT64_TO_FRAC(-1);

	PPG_FRAC max_xj = ZERO;
	PPG_FRAC max_dxj = ZERO;

	// Calculate A_norm2
	for (size_t j = 0; j < N; ++j) {
		A_norm2[j] = ZERO;
		for (size_t i = 0; i < M; ++i) {
			A_norm2[j] = ADD(A_norm2[j], MULT(A[i * N + j], A[i * N + j]));
		}
	}

	// Init x_hat
	for (size_t j = 0; j < N; ++j) {
		x_hat[j] = ZERO;
	}

	// Init max_xj: could have just set to 0.5, but if we change to random
	// init or something, then this can be useful.
	for (size_t j = 0; j < N; ++j) {
		max_xj = MAX(FABS(x_hat[j]), max_xj);
	}

	// Set up residual
	for (size_t i = 0; i < M; ++i) {
		r[i] = y[i];
		for (size_t j = 0; j < N; ++j) {
			r[i] = SUB(r[i], MULT(A[i * N + j], x_hat[j]));
		}
	}

	do {
		max_dxj = ZERO;
		for (size_t j = 0; j < N; ++j) {
			// A_norm2[j] needs to be non-zero
			if (A_norm2[j] == ZERO) continue;
			// keep old value
			PPG_FRAC x_j0 = x_hat[j];
			// Set up x[j] update
			PPG_FRAC rho_j = ZERO;
			for (size_t i = 0; i < M; ++i) {
				// remove x[j] contribution from residual
				r[i] = ADD(r[i], MULT(A[i * N + j], x_hat[j]));
				rho_j = ADD(rho_j, MULT(r[i], A[i * N + j]));
			}
			PPG_FRAC sign = GREATER_THAN(rho_j, ZERO) ? ONE : MINUS_ONE;
			x_hat[j] = DIV(MULT(sign, MAX(SUB(FABS(rho_j), lambda),
																		ZERO)),
										 A_norm2[j]);
			for (size_t i = 0; i < M; ++i) {
				// add back x[j] contribution to residual, using new value
				r[i] = SUB(r[i], MULT(A[i * N + j], x_hat[j]));
			}
			max_dxj = MAX(FABS(SUB(x_hat[j], x_j0)), max_dxj);
			max_xj = MAX(FABS(x_hat[j]), max_xj);
		}
	} while (GREATER_THAN(DIV(max_dxj, max_xj), tol));

	free(r);
	free(A_norm2);
}

// PPG algorithm helper functions
void setup_psiT(PPG_FRAC *psiT, size_t N) {
	PPG_FRAC factor = FLOAT64_TO_FRAC(1.0 / sqrt(N));
	for (size_t n = 0; n < N; ++n) {
		psiT[n * N] = factor;
	}

	factor = MULT(factor, FLOAT64_TO_FRAC(sqrt(2)));
	for (size_t n = 0; n < N; ++n) {
		for (size_t k = 1; k < N; ++k) {
			psiT[n * N + k] = MULT(factor,
														 FLOAT64_TO_FRAC(cos((M_PI / N) * (n + 0.5) * k)));
		}
	}
}

void half_flip(bool *flags, size_t N, size_t M, bool *flipped) {
	size_t m = 0;
	for (size_t i = 0; i < N; ++i) {
		flipped[i] = false;
	}

	for (size_t i = 0; i < N; ++i) {
		if (flags[i]) {
			if (m++ < M/2) {
				// should be moved forward
				// but only if possible...
				if (i < N/2) {
					flipped[i+N/2] = true;
				}
			} else {
				// should be moved backward
				// but only if possible...
				if (i >= N/2) {
					flipped[i-N/2] = true;
				}
			}
		}
		if (m == M) break;
	}
}

void get_selected_rows(PPG_FRAC *mat, bool *flags, 
											size_t N, size_t m, PPG_FRAC *out) {
	size_t row = 0;
	for (size_t i = 0; i < N; ++i) {
		if (flags[i]) {
			memcpy(out + row * N, mat + i * N, N * sizeof(PPG_FRAC));
			if (row++ == m) break;
		}
	}
}

void ppg(PPG_FRAC *Y, bool *phi_flags, PPG_Params params,
				PPG_FRAC *tr, PPG_FRAC *Xr) {
	PPG_FRAC *y = malloc(sizeof(PPG_FRAC) * params.N_window);
	PPG_FRAC *s = malloc(sizeof(PPG_FRAC) * params.N_window);
	PPG_FRAC *xr = malloc(sizeof(PPG_FRAC) * params.N_window);

	// IDCT matrix
	PPG_FRAC *psiT 
			= malloc(sizeof(PPG_FRAC) * params.N_window * params.N_window);
	setup_psiT(psiT, params.N_window);

	// A = phi * psi.T
	PPG_FRAC *A = malloc(sizeof(PPG_FRAC) * params.M * params.N_window);
	get_selected_rows(psiT, phi_flags, params.N_window, params.M, A);

	// Flipped A for intermediate windows.
	bool *phi_flags_flip = malloc(sizeof(bool) * params.N_window);
	PPG_FRAC *A_flip = malloc(sizeof(PPG_FRAC) * params.M * params.N_window);
	half_flip(phi_flags, params.N_window, params.M, phi_flags_flip);
	get_selected_rows(psiT, phi_flags_flip, params.N_window, params.M,
										A_flip);

	size_t N0 = params.N0;
	size_t M0 = N0 / params.CF;
	size_t M = params.M;
	size_t N_window = params.N_window;

	// Setup reconstruction times first, as we know those (uniform
	// sampling)
	for (size_t i = 0; i < N0; ++i) {
		tr[i] = MULT(INT64_TO_FRAC(i), params.T0);
	}

	// Standardize scale
#if defined(PPG_USE_FXP64) || defined(PPG_USE_FXP32) \
		|| defined(PPG_USE_FXP16) || defined(PPG_USE_FXP8)
	standardize(Y, M0, 1 << FxP_INT_LEN);
#endif

	// First quarter-window
	memcpy(y, Y, sizeof(PPG_FRAC) * M);
	cd_lasso(y, A, M, N_window, params.LASSO_lambda,
					params.CD_diff_thresh, s);
	dot(psiT, s, N_window, N_window, 1, xr);
	memcpy(Xr, xr, sizeof(PPG_FRAC) * N_window / 4);
	
	// All the windows
	for (size_t t_x = 0, t_y = 0;
			 t_y < M0 - M + 1;
			 t_x += N_window / 2, t_y += M / 2) {
		memcpy(y, Y + t_y, sizeof(PPG_FRAC) * M);
		if (t_y % M == 0) {
			cd_lasso(y, A, M, N_window,
								params.LASSO_lambda, params.CD_diff_thresh,
								s);
		} else {
			cd_lasso(y, A_flip, M, N_window, params.LASSO_lambda,
							params.CD_diff_thresh, s);
		}
		dot(psiT, s, N_window, N_window, 1, xr);
		memcpy(Xr + t_x + N_window / 4, xr + N_window / 4,
					 sizeof(PPG_FRAC) * N_window / 2);

		// TODO: add max_f in setup somewhere, instead of hardcoding like
		// this.
		double bpm_window = 60.0 * get_freq(Xr + t_x + N_window / 4,
																				tr + t_x + N_window / 4,
																				N_window / 2,
																				4.0);
		printf("BPM for (%.2f, %.2f): %.2f\n", FRAC_TO_FLOAT64(tr[t_x + N_window/4]),
																					 FRAC_TO_FLOAT64(tr[t_x + 3*N_window/4]),
																					 bpm_window);
	}

	// Last quarter window
	memcpy(y, Y + M0 - M, sizeof(double) * M);
	if (M0 % M == 0) {
		cd_lasso(y, A, M, N_window, params.LASSO_lambda,
						params.CD_diff_thresh, s);
	} else {
		cd_lasso(y, A_flip, M, N_window, params.LASSO_lambda,
						params.CD_diff_thresh, s);
	}
	dot(psiT, s, N_window, N_window, 1, xr);
	memcpy(Xr + N0 - N_window/4, xr + 3 * N_window/4,
				sizeof(double) * N_window / 4);

	free(phi_flags_flip);
	free(A_flip);
	free(A);
	free(psiT);
	free(y);
	free(xr);
	free(s);
}

// see py/main.py for more intuitive code, and comments
double get_freq(PPG_FRAC *X, PPG_FRAC *t, size_t N, double max_f) {
	double t_min = 1.0 / max_f;
	double mean_X = 0;
	for (size_t i = 0; i < N; ++i) {
		mean_X += (1.0 / N) * FRAC_TO_FLOAT64(X[i]);
	}
	
	double rms = 0.0;
	for (size_t i = 0; i < N; ++i) {
		rms += (1.0 / N) * (FRAC_TO_FLOAT64(X[i]) - mean_X) 
										 * (FRAC_TO_FLOAT64(X[i]) - mean_X);
	}
	rms = sqrt(rms);
	
	double xprod;
	double xsum;
	size_t last_cross_positive = 0;
	size_t last_cross_negative = 0;
	size_t nc = 0;
	for (size_t i = 0; i < N-1; ++i) {
		// Looking for positive crossing
		xprod = (FRAC_TO_FLOAT64(X[i]) - mean_X - rms) 
					* (FRAC_TO_FLOAT64(X[i+1]) - mean_X - rms);
		if (xprod <= 0) {
			xsum = fabs(FRAC_TO_FLOAT64(X[i]) - mean_X - rms) 
						+ fabs(FRAC_TO_FLOAT64(X[i+1]) - mean_X - rms);
			if ((xsum != 0.0) 
					&& ((FRAC_TO_FLOAT64(t[i]) 
								- FRAC_TO_FLOAT64(t[last_cross_positive])) 
							>= t_min)) {
				nc++;
				last_cross_positive = i;
			}
			continue; // should probably move on to the next sample now
		}

		// Negative crossings now.
		xprod = (FRAC_TO_FLOAT64(X[i]) - mean_X + rms) 
					* (FRAC_TO_FLOAT64(X[i+1]) - mean_X + rms);
		if (xprod <= 0) {
			xsum = fabs(FRAC_TO_FLOAT64(X[i]) - mean_X + rms) 
						+ fabs(FRAC_TO_FLOAT64(X[i+1]) - mean_X + rms);
			if ((xsum != 0.0) 
					&& ((FRAC_TO_FLOAT64(t[i])
								- FRAC_TO_FLOAT64(t[last_cross_negative]))
							>= t_min)) {
				nc++;
				last_cross_negative = i;
			}
		}
	}
	return ((double) nc) 
					/ (2 * (FRAC_TO_FLOAT64(t[N-1]) - FRAC_TO_FLOAT64(t[0])));
}

double corrcoef(PPG_FRAC *a, PPG_FRAC *b, size_t N) {
	double mean_a, mean_b;
	mean_a = mean_b = 0.0;
	for (size_t i = 0; i < N; ++i) {
		mean_a += (1.0 / N) * FRAC_TO_FLOAT64(a[i]);
		mean_b += (1.0 / N) * FRAC_TO_FLOAT64(b[i]);
	}

	double norm_a = 0.0, norm_b = 0.0;
	for (size_t i = 0; i < N; ++i) {
		norm_a += pow(FRAC_TO_FLOAT64(a[i]) - mean_a, 2);
		norm_b += pow(FRAC_TO_FLOAT64(b[i]) - mean_b, 2);
	}
	norm_a = sqrt(norm_a);
	norm_b = sqrt(norm_b);

	// Don't let (1/std) blow up; better to just let it die.
	if (norm_a < 1.0e-4) norm_a = 1.0e10;
	if (norm_b < 1.0e-4) norm_b = 1.0e10;

	double corr = 0.0;
	for (size_t i = 0; i < N; ++i) {
		corr += ((FRAC_TO_FLOAT64(a[i]) - mean_a) 
						* (FRAC_TO_FLOAT64(b[i]) - mean_b));
	}
	corr /= (norm_a * norm_b);
	return corr;
}

