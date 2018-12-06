#include "ppg_algorithm.h"
#include <math.h>
#include <string.h>
#include <stdio.h>

// General helper functions
static double max(double a, double b) {
	return (a > b ? a : b) ;
}

static void vec_sub_update(double *vec, double *b, size_t N) {
	for (size_t i = 0; i < N; ++i) {
		vec[i] -= b[i];
	}
}

////
// This portion is for the code I was writing to use the DCT definition
// directly instead of storing the cosine values in a matrix. This could
// be useful later, so I am not removing it yet.
//
// TODO: figure out how to not lose this code using the power of GIT!!!
//
// DCT stuff
// O(N^2) algorithms: just naive loops
void dct_1D(double *X, size_t N, double *X_DCT) {
	memset(X_DCT, 0, sizeof(double) * N);
	const double factor_1 = 2.0 / sqrt(N);
	const double factor_0 = factor_1 / sqrt(2);
	for (size_t k = 0; k < N; ++k) {
		for (size_t n = 0; n < N; ++n) {
			X_DCT[k] += cos((M_PI / N) * (n + 0.5) * k) * X[n];
		}
	}
	X_DCT[0] *= factor_0;
	for (size_t k = 1; k < N; ++k) {
		X_DCT[k] *= factor_1;
	}
}

void idct_1D(double *X_DCT, size_t N, double *X) {
	memset(X, 0, sizeof(double) * N);
	const double factor_1 = 2.0 / sqrt(N);
	const double factor_0 = factor_1 / sqrt(2);
	for (size_t n = 0; n < N; ++n) {
		X[n] += X_DCT[0] * factor_0;
	}

	for (size_t k = 1; k < N; ++k) {
		for (size_t n = 0; n < N; ++n) {
			X[n] += cos((M_PI / N) * (n + 0.5) * k) * X_DCT[k] * factor_1;
		}
	}
}

void random_sample_idct_1D(double *X_DCT, size_t N,
													bool *phi_flags, size_t M, double *Y) {
	memset(Y, 0, sizeof(double) * M);
	const double factor_1 = 2.0 / sqrt(N);
	const double factor_0 = factor_1 / sqrt(2);
	for (size_t m = 0; m < M; ++m) {
		Y[m] = X_DCT[0] * factor_0;
	}

	size_t m = 0;
	for (size_t n = 0; n < N; ++n) {
		if (phi_flags[n]) {
			for (size_t k = 1; k < N; ++k) {
				Y[m] += cos((M_PI / N) * (n + 0.5) * k) * X_DCT[k] * factor_1;
			}
			m++;
		}
	}
}

// LASSO solver helper functions
static void calc_Anorm2_randsample_idct_1D(bool *phi_flags, size_t N,
																					double *A_norm2) {
	memset(A_norm2, 0, sizeof(double) * N);
	const double factor_1 = 2.0 / sqrt(N);
	const double factor_0 = factor_1 / sqrt(2);
	for (size_t n = 0; n < N; ++n) {
		if (phi_flags[n]) {
			A_norm2[0] += factor_0 * factor_0;
			for (size_t k = 1; k < N; ++k) {
				A_norm2[k] += pow(cos((M_PI / N) * (n + 0.5) * k) * factor_1, 2);
			}
		}
	}
}

void cd_lasso_randsample_idct_1D(double *y, size_t M,
																bool *phi_flags, size_t N,
																double lambda, double tol,
																double *s_hat) {
	// Y is M-dim, phi_flags controls which DCT components to use, x_hat
	// is Nx1 
	// We don't need an explicit (A) matrix, as we know it's going to
	// randomly sample from the IDCT. We can call the functions which can
	// calculate the products we want from phi_flags and the definition of
	// the DCT.
	double *A_norm2 = malloc(sizeof(double) * N);
	calc_Anorm2_randsample_idct_1D(phi_flags, N, A_norm2);
	
	// can initialize with anything I think...
	for (size_t k = 0; k < N; ++k) {
		s_hat[k] = 0.5;
	}

	double *y_hat = malloc(sizeof(double) * M);
	random_sample_idct_1D(s_hat, N, phi_flags, M, y_hat);

	double *r = malloc(sizeof(double) * M);
	memcpy(r, y, sizeof(double) * M);
	vec_sub_update(r, y_hat, M);

	const double factor_1 = 2.0 / sqrt(N);
	const double factor_0 = factor_1 / sqrt(2);
	
	double max_sk = 0.5; // depends on initialization
	printf("CD!\n");
	while (true) {
		double max_dsk = 0;
		for (size_t k = 0; k < N; ++k) {
			if (A_norm2[k] == 0.0)
				continue;

			double factor = (k == 0 ? factor_0 : factor_1);
			double s_k0 = s_hat[k];

			// r += A[:,k] .* s[k]
			double rho_k = 0;
			size_t m = 0;
			for (size_t n = 0; n < N; ++n) {
				if (phi_flags[n]) {
					r[m] += s_hat[k] * cos((M_PI / N) * (n + 0.5) * k) * factor;
					rho_k += r[m] * cos((M_PI / N) * (n + 0.5) * k) * factor;
					m++;
				}
			}

			double sign = (rho_k > 0 ? 1.0 : -1.0);
			s_hat[k] = (sign * max(fabs(rho_k) - lambda, 0)) / A_norm2[k];

			double dsk = fabs(s_hat[k] - s_k0);
			max_dsk = max(dsk, max_dsk);
			max_sk = max(fabs(s_hat[k]), max_sk);

			// Adding back the residual contribution from s_hat[k]
			m = 0;
			for (size_t n = 0; n < N; ++n) {
				if (phi_flags[n]) {
					r[m++] -= s_hat[k] * cos((M_PI / N) * (n + 0.5) * k) * factor;
				}
			}
		}
		printf("%lf %lf\n", max_sk, max_dsk);
		if ((max_dsk / max_sk) < tol) break;
	}

	free(r);
	free(y_hat);
	free(A_norm2);
}

// END OF DCT stuff
////

////////////////////////////////////////////////////////////////////////
// THIS IS THE REAL SHIT. THE CODE WHICH IS ACTUALLY RUNNING.
// Matrix math
void dot(double *A, double *B, size_t M, size_t N, size_t P, double *C) {
	// C = AB, everything is row-major
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < P; ++j) {
			C[i * P + j] = 0;
			for (int k = 0; k < N; ++k) {
				C[i * P + j] += A[i * N + k] * B[k * P + j];
			}
		}
	}
}

// LASSO solvers
void cd_lasso(double *y, double *A, size_t M, size_t N,
				double lambda, double tol, double *x_hat) {
	// y M-dim, A is MxN, x_hat is Nx1
	double *A_norm2 = malloc(sizeof(double) * N);
	memset(A_norm2, 0, sizeof(double) * N);
	for (int i = 0; i < M; ++i) {
		for (int j = 0; j < N; ++j) {
			A_norm2[j] += pow(A[i * N + j], 2);
		}
	}

	for (int i = 0; i < N; ++i) {
		x_hat[i] = 0.5;
	}

	double *y_hat = malloc(sizeof(double) * M);
	// y_hat = A * x_hat
	dot(A, x_hat, M, N, 1, y_hat);

	double *r = malloc(sizeof(double) * M);
	memcpy(r, y, sizeof(double) * M);
	vec_sub_update(r, y_hat, M);

	double max_xi = 0;
	for (int i = 0; i < N; ++i) {
		if (fabs(x_hat[i]) > max_xi) {
			max_xi = fabs(x_hat[i]);
		}
	}

	while (true) {
		double max_dxi = 0;
		for (int i = 0; i < N; ++i) {
			if (A_norm2[i] == 0.0) {
				continue;
			}

			double x_i0 = x_hat[i];

			// r += A[:,i] .* x[i]
			double rho_i = 0;
			for (int j = 0; j < M; ++j) {
				r[j] += A[j * N + i] * x_hat[i];
				rho_i += r[j] * A[j * N + i];
			}

			double sign = (rho_i > 0 ? 1.0 : -1.0);
			x_hat[i] = (sign * max(fabs(rho_i) - lambda, 0)) / A_norm2[i];

			double dxi = fabs(x_hat[i] - x_i0);
			max_dxi = (dxi > max_dxi ? dxi : max_dxi);

			for (int j = 0; j < M; ++j) {
				r[j] -= x_hat[i] * A[j * N + i];
			}

			max_xi = max(max_xi, fabs(x_hat[i]));
		}
		if ((max_dxi / max_xi) < tol) {
			break;
		}
	}

	free(r);
	free(y_hat);
	free(A_norm2);
}

// PPG algorithm helper functions
static void setup_psiT(double *psiT, size_t N) {
	double factor = 2.0 / sqrt(2 * N);
	for (size_t n = 0; n < N; ++n) {
		psiT[n * N] = factor;
	}

	factor *= sqrt(2);
	for (size_t n = 0; n < N; ++n) {
		for (size_t k = 0; k < N; ++k) {
			psiT[n * N + k] = factor * cos((M_PI / N) * (n + 0.5) * k);
		}
	}
}

static void half_flip(bool *flags, size_t N, size_t M, bool *flipped) {
	size_t m = 0;
	for (size_t i = 0; i < N; ++i) {
		if (flags[i]) {
			if (m++ < M/2) {
				// should be moved forward
				flipped[i+N/2] = true;
			} else {
				// should be moved backward
				flipped[i-N/2] = true;
			}
		}
		if (m == M) break;
	}
}

static void get_selected_rows(double *mat, bool *flags, 
															size_t N, size_t m, double *out) {
	size_t row = 0;
	for (int i = 0; i < N; ++i) {
		if (flags[i]) {
			memcpy(out + row * N, mat + i * N, N * sizeof(double));
			if (row++ == m) break;
		}
	}
}

void ppg(double *Y, bool *phi_flags, PPG_Params params, double *Xr) {
	/*
	 * t0		:	sample instants (in s) 
	 * X0		:	original samples (which we are trying to reconstruct)
	 * N0		:	number of samples we have
	 * phi_flags:	boolean flags to denote which observations we are selecting,
	 *				this should have been randomly generated previously
	 * params	:	simulation/experiment parameters
	 * Xr		:	where we should store the reconstructed signal (output)
	 */
	
	double *y = malloc(sizeof(double) * params.N_window);
	double *s = malloc(sizeof(double) * params.N_window);
	double *xr = malloc(sizeof(double) * params.N_window);

	// IDCT matrix
	double *psiT 
			= malloc(sizeof(double) * params.N_window * params.N_window);
	setup_psiT(psiT, params.N_window);

	// A = phi * psi.T
	double *A = malloc(sizeof(double) * params.M * params.N_window);
	get_selected_rows(psiT, phi_flags, params.N_window, params.M, A);

	// Flipped A for intermediate windows.
	// TODO: explain properly why this is needed.
	bool *phi_flags_flip = malloc(sizeof(bool) * params.N_window);
	double *A_flip = malloc(sizeof(double) * params.M * params.N_window);
	half_flip(phi_flags, params.N_window, params.M, phi_flags_flip);
	get_selected_rows(psiT, phi_flags_flip, params.N_window, params.M,
										A_flip);

	size_t N0 = params.N0;
	size_t M0 = N0 / params.CF;
	size_t M = params.M;
	size_t N_window = params.N_window;

	// First quarter-window
	memcpy(y, Y, sizeof(double) * M);
	cd_lasso(y, A, M, N_window, params.LASSO_lambda,
					params.CD_diff_thresh, s);
	dot(psiT, s, N_window, N_window, 1, xr);
	memcpy(Xr, xr, sizeof(double) * N_window / 4);
	
	// All the windows
	for (size_t t_x = 0, t_y = 0;
			 t_y < M0 - M + 1;
			 t_x += N_window / 2, t_y += M / 2) {
		memcpy(y, Y + t_y, sizeof(double) * M);
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
					 sizeof(double) * N_window / 2);
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