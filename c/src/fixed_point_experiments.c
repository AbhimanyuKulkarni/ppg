/*
 * fixed_point_experiments.c
 *
 * Some basic fiddling with fixed point number representations.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "fixed_point.h"

static double max_double(double a, double b) {
	return a > b ? a : b;
}

static FxP64 max_FxP64(FxP64 a, FxP64 b) {
	// signed ints should work the same
	return a > b ? a : b;
}

static FxP64 abs_FxP64(FxP64 x) {
	// signed ints should work the same?
	return x < 0 ? -x : x;
}

static bool less_than_FxP64(FxP64 a, FxP64 b) {
	return sub_FxP64(a, b) < 0;
}

static bool greater_than_FxP64(FxP64 a, FxP64 b) {
	return sub_FxP64(a, b) > 0;
}

void test_FxP64_conversion(double *nums, size_t N) {
	printf("=========================================================\n");
	printf("Testing conversion:\n");
	for (size_t i = 0; i < N; ++i) {
		printf("=================================\n");
		union double_bitview in;
		in.dval = nums[i];

		FxP64 fxp_i = double_to_FxP64(in.dval);
		union double_bitview out;
		out.dval = FxP64_to_double(fxp_i);
		if (out.dval == in.dval) {
			printf("(OK) ");
		} else {
			printf("(Failed) ");
		}
		printf("Double %.10lf (%016lx) -> %016lx; ",
						in.dval, in.uintval, fxp_i);
		printf("And %016lx -> %.10lf (%016lx)\n",
						fxp_i, out.dval, out.uintval);
	}
	printf("=========================================================\n");
}

void test_FxP64_arithmetic(double *nums, size_t N, double tol) {
	printf("=========================================================\n");
	printf("Testing arithmetic:\n");
	double *sums = malloc(sizeof(double) * N * N);
	double *diffs = malloc(sizeof(double) * N * N);
	double *prods = malloc(sizeof(double) * N * N);
	double *divs = malloc(sizeof(double) * N * N);
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			sums[i*N+j] = nums[i] + nums[j];
			diffs[i*N+j] = nums[i] - nums[j];
			prods[i*N+j] = nums[i] * nums[j];
			if (nums[j] != 0.0) {
				divs[i*N+j] = nums[i] / nums[j];
			} else {
				divs[i*N+j] = 0.0;
			}
		}
	}

	for (size_t i = 0; i < N; ++i) {
		printf("=================================\n");
		FxP64 fxp_i = double_to_FxP64(nums[i]);
		for (size_t j = 0; j < N; ++j) {
			FxP64 fxp_j = double_to_FxP64(nums[j]);
			double sum = FxP64_to_double(add_FxP64(fxp_i, fxp_j));
			double diff = FxP64_to_double(sub_FxP64(fxp_i, fxp_j));
			double prod = FxP64_to_double(mult_FxP64(fxp_i, fxp_j));
			double div = 0.0;
			if (fxp_j != 0) {
				div = FxP64_to_double(div_FxP64(fxp_i, fxp_j));
			}
			
			if (fabs(sum - sums[i*N+j]) > tol) {
				printf("Failed");
			} else {
				printf("OK");
			}
			printf("\t%lf + %lf = %lf, expected %lf\n", nums[i], nums[j],
																									sum, sums[i*N+j]);

			if (fabs(diff - diffs[i*N+j]) > tol) {
				printf("Failed");
			} else {
				printf("OK");
			}
			printf("\t%lf - %lf = %lf, expected %lf\n", nums[i], nums[j],
																									diff, diffs[i*N+j]);

			if (fabs(prod - prods[i*N+j]) > tol) {
				printf("Failed");
			} else {
				printf("OK");
			}
			printf("\t%lf * %lf = %lf, expected %lf\n", nums[i], nums[j],
																									prod, prods[i*N+j]);

			if (fabs(div - divs[i*N+j]) > tol) {
				printf("Failed");
			} else {
				printf("OK");
			}
			printf("\t%lf / %lf = %lf, expected %lf\n", nums[i], nums[j],
																									div, divs[i*N+j]);
		}
	}

	free(sums);
	free(diffs);
	free(prods);
	free(divs);
	printf("=========================================================\n");
}

void test_FxP64_dot_products(size_t p_max) {
	printf("=========================================================\n");
	printf("Testing dot Products:\n");
	for (size_t p = 0; p < p_max; ++p) {
		double *vec_a = malloc(sizeof(double) * (1 << p));
		double *vec_b = malloc(sizeof(double) * (1 << p));

		double num;
		for (size_t i = 0; i < (1 << p); ++i) {
			num = rand() % (1 << (FxP64_INT_LEN - 1));
			num /= (1 << (FxP64_INT_LEN - 1)); // in [0, 1]
			num /= sqrt(1 << p);
			vec_a[i] = num;
		}
		for (size_t i = 0; i < (1 << p); ++i) {
			num = rand() % (1 << (FxP64_INT_LEN - 1));
			num /= (1 << (FxP64_INT_LEN - 1)); // in [0, 1]
			num /= sqrt(1 << p);
			vec_b[i] = num;
		}

		double dot = 0.0;
		FxP64 dot_FxP64 = double_to_FxP64(0.0);
		for (size_t i = 0; i < (1 << p); ++i) {
			dot += vec_a[i] * vec_b[i];
			dot_FxP64 = add_FxP64(dot_FxP64, mult_FxP64(double_to_FxP64(vec_a[i]),
																									double_to_FxP64(vec_b[i])));
		}

		printf("p = %lu; double dot: %.6lf, FxP64 dot: %.6lf\n",
					p, dot, FxP64_to_double(dot_FxP64));

		free(vec_a);
		free(vec_b);
	}
	printf("=========================================================\n");
}

void cd_lasso_double(double *y, double *A, double *x,
										size_t N, size_t D,
										double lambda, double tol) {
	double *A_norm2 = malloc(sizeof(double) * D);
	for (size_t j = 0; j < D; ++j) {
		A_norm2[j] = 0.0;
		for (size_t i = 0; i < N; ++i) {
			A_norm2[j] += A[i * D + j] * A[i * D + j];
		}
	}

	double max_xj = 0;
	for (size_t j = 0; j < D; ++j) {
		max_xj = max_double(fabs(x[j]), max_xj);
	}

	double *r = malloc(sizeof(double) * N);
	for (size_t i = 0; i < N; ++i) {
		r[i] = y[i];
		for (size_t j = 0; j < D; ++j) {
			r[i] -= A[i * D + j] * x[j];
		}
	}

	double max_dxj = 0.0;
	do {
		max_dxj = 0.0;
		for (size_t j = 0; j < D; ++j) {
			if (A_norm2[j] == 0.0) continue;
			double x_j0 = x[j];
			double rho_j = 0.0;
			for (size_t i = 0; i < N; ++i) {
				r[i] += A[i * D + j] * x[j];
				rho_j += r[i] * A[i * D + j];
			}
			double sign = (rho_j > 0 ? 1.0 : -1.0);
			x[j] = (sign * max_double(fabs(rho_j) - lambda, 0)) / A_norm2[j];
			for (size_t i = 0; i < N; ++i) {
				r[i] -= A[i * D + j] * x[j];
			}
			max_dxj = max_double(fabs(x[j] - x_j0), max_dxj);
			max_xj = max_double(fabs(x[j]), max_xj);
		}
	} while ((max_dxj / max_xj) < tol);
	free(r);
	free(A_norm2);
}

void cd_lasso_FxP64(FxP64 *y, FxP64 *A, FxP64 *x, size_t N, size_t D,
										FxP64 lambda, FxP64 tol) {
	FxP64 *A_norm2 = malloc(sizeof(FxP64) * D);
	const FxP64 ZERO = int64_to_FxP64((int64_t) 0);
	const FxP64 ONE = int64_to_FxP64((int64_t) 1);
	const FxP64 MINUS_ONE = int64_to_FxP64((int64_t) -1);

	for (size_t j = 0; j < D; ++j) {
		A_norm2[j] = ZERO;
		for (size_t i = 0; i < N; ++i) {
			A_norm2[j] = add_FxP64(A_norm2[j], mult_FxP64(A[i * D + j],
																										A[i * D + j]));
		}
	}

	FxP64 max_xj = ZERO;
	for (size_t j = 0; j < D; ++j) {
		max_xj = max_FxP64(abs_FxP64(x[j]), max_xj);
	}

	FxP64 *r = malloc(sizeof(FxP64) * N);
	for (size_t i = 0; i < N; ++i) {
		r[i] = y[i];
		for (size_t j = 0; j < D; ++j) {
			r[i] = sub_FxP64(r[i], mult_FxP64(A[i * D + j], x[j]));
		}
	}

	FxP64 max_dxj = ZERO;
	do {
		max_dxj = ZERO;
		for (size_t j = 0; j < D; ++j) {
			if (A_norm2[j] == ZERO) continue;
			FxP64 x_j0 = x[j];
			FxP64 rho_j = ZERO;
			for (size_t i = 0; i < N; ++i) {
				r[i] = add_FxP64(r[i], mult_FxP64(A[i * D + j], x[j]));
				rho_j = add_FxP64(rho_j, mult_FxP64(r[i], A[i * D + j]));
			}
			FxP64 sign = greater_than_FxP64(rho_j, ZERO) ? ONE : MINUS_ONE;
			x[j] = div_FxP64(mult_FxP64(sign,
																	max_FxP64(sub_FxP64(abs_FxP64(rho_j),
																											lambda),
																						ZERO)),
											A_norm2[j]);
			for (size_t i = 0; i < N; ++i) {
				r[i] = sub_FxP64(r[i], mult_FxP64(A[i * D + j], x[j]));
			}
			max_dxj = max_FxP64(abs_FxP64(sub_FxP64(x[j], x_j0)), max_dxj);
			max_xj = max_FxP64(abs_FxP64(x[j]), max_xj);
		}
	} while (less_than_FxP64(div_FxP64(max_dxj, max_xj), tol));

	free(r);
	free(A_norm2);
}

void test_lasso(size_t N) {
	printf("=========================================================\n");
	printf("Testing with LASSO:\n");
	double *A_double = malloc(sizeof(double) * N * 2);
	double *y_double = malloc(sizeof(double) * N);
	double *x_double = malloc(sizeof(double) * 2);
	x_double[0] = 0.0;
	x_double[1] = 1.0;
	printf("True parameter: (%lf, %lf)\n", x_double[0], x_double[1]);

	for (size_t i = 0; i < N; ++i) {
		A_double[i*2] = 2 * (((double) rand()) / RAND_MAX - 0.5);
		A_double[i*2 + 1] = 2 * (((double) rand()) / RAND_MAX - 0.5);
		y_double[i] = A_double[i*2] * x_double[0]
									+ A_double[i*2 + 1] * x_double[1];
	}
	x_double[0] = 0.5;
	x_double[1] = 0.5;
	cd_lasso_double(y_double, A_double, x_double, N, 2, 0.1, 0.0001);
	printf("Estimated parameter (using double): (%lf, %lf)\n",
					x_double[0], x_double[1]);

	FxP64 *A, *y, *x;
	A = malloc(sizeof(FxP64) * N * 2);
	y = malloc(sizeof(FxP64) * N);
	x = malloc(sizeof(FxP64) * 2);
	x[0] = double_to_FxP64(0.5);
	x[1] = double_to_FxP64(0.5);
	for(size_t i = 0; i < N; ++i) {
		A[i*2] = double_to_FxP64(A_double[i*2]);
		A[i*2+1] = double_to_FxP64(A_double[i*2+1]);
		y[i] = double_to_FxP64(y_double[i]);
	}
	cd_lasso_FxP64(y, A, x, N, 2,
								double_to_FxP64(0.1), double_to_FxP64(0.0001));
	printf("Estimated parameter (using FxP64): (%lf, %lf)\n",
													FxP64_to_double(x[0]), FxP64_to_double(x[1]));

	free(A);
	free(y);
	free(x);
	free(y_double);
	free(x_double);
	free(A_double);
	printf("=========================================================\n");
}

int main(int argc, char *argv[]) {
	const size_t N = 32;
	double nums[32] = {0.0, -0.0, 1.0, -1.0,
									2.0, 8.0, 32767.0, -32767.0,
									0.5, 0.25, 0.125, 0.0625,
									1.0e-1, 1.0e-4, 1.0e-8, 1.0e-16,
									1.0e1, 1.0e2, 1.0e4, 1.0e8,
									1.0/3, M_PI, -1.0/3, -1.25,
									-1.0e-1, -1.0e-4, -1.0e-8, -1.0e-16,
									-1.0e1, -1.0e2, -1.0e4, -1.0e8
								 };
	test_FxP64_conversion(nums, N);
	const double tol = 1.0e-6;
	test_FxP64_arithmetic(nums, N, tol);
	test_FxP64_dot_products(20);
	srand(0);
	test_lasso(30);
	return 0;
}
