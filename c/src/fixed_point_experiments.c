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

#include "fxp_ops.h"

static double max_double(double a, double b) {
	return a > b ? a : b;
}

void test_FxP_conversion(double *nums, size_t N) {
	printf("=========================================================\n");
	printf("Testing conversion:\n");
	for (size_t i = 0; i < N; ++i) {
		printf("=================================\n");
		union double_bitview in;
		in.dval = nums[i];

		FxP fxp_i = double_to_FxP(in.dval);
		union double_bitview out;
		out.dval = FxP_to_double(fxp_i);
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

void test_FxP_arithmetic(double *nums, size_t N, double tol) {
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
		FxP fxp_i = double_to_FxP(nums[i]);
		for (size_t j = 0; j < N; ++j) {
			FxP fxp_j = double_to_FxP(nums[j]);
			double sum = FxP_to_double(add_FxP(fxp_i, fxp_j));
			double diff = FxP_to_double(sub_FxP(fxp_i, fxp_j));
			double prod = FxP_to_double(mult_FxP(fxp_i, fxp_j));
			double div = 0.0;
			if (fxp_j != 0) {
				div = FxP_to_double(div_FxP(fxp_i, fxp_j));
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

void test_FxP_dot_products(size_t p_max) {
	printf("=========================================================\n");
	printf("Testing dot Products:\n");
	for (size_t p = 0; p < p_max; ++p) {
		double *vec_a = malloc(sizeof(double) * (1 << p));
		double *vec_b = malloc(sizeof(double) * (1 << p));

		double num;
		for (size_t i = 0; i < (1 << p); ++i) {
			num = rand() % (1 << (FxP_INT_LEN - 1));
			num /= (1 << (FxP_INT_LEN - 1)); // in [0, 1]
			num /= sqrt(1 << p);
			vec_a[i] = num;
		}
		for (size_t i = 0; i < (1 << p); ++i) {
			num = rand() % (1 << (FxP_INT_LEN - 1));
			num /= (1 << (FxP_INT_LEN - 1)); // in [0, 1]
			num /= sqrt(1 << p);
			vec_b[i] = num;
		}

		double dot = 0.0;
		FxP dot_FxP = double_to_FxP(0.0);
		for (size_t i = 0; i < (1 << p); ++i) {
			dot += vec_a[i] * vec_b[i];
			dot_FxP = add_FxP(dot_FxP, mult_FxP(double_to_FxP(vec_a[i]),
																					double_to_FxP(vec_b[i])));
		}

		printf("p = %lu; double dot: %.6lf, FxP dot: %.6lf\n",
					p, dot, FxP_to_double(dot_FxP));

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
	} while ((max_dxj / max_xj) > tol);
	free(r);
	free(A_norm2);
}

void cd_lasso_FxP(FxP *y, FxP *A, FxP *x, size_t N, size_t D,
										FxP lambda, FxP tol) {
	FxP *A_norm2 = malloc(sizeof(FxP) * D);
	const FxP ZERO = int64_to_FxP((int64_t) 0);
	const FxP ONE = int64_to_FxP((int64_t) 1);
	const FxP MINUS_ONE = int64_to_FxP((int64_t) -1);

	for (size_t j = 0; j < D; ++j) {
		A_norm2[j] = ZERO;
		for (size_t i = 0; i < N; ++i) {
			A_norm2[j] = add_FxP(A_norm2[j], mult_FxP(A[i * D + j],
																										A[i * D + j]));
		}
	}

	FxP max_xj = ZERO;
	for (size_t j = 0; j < D; ++j) {
		max_xj = max_FxP(abs_FxP(x[j]), max_xj);
	}

	FxP *r = malloc(sizeof(FxP) * N);
	for (size_t i = 0; i < N; ++i) {
		r[i] = y[i];
		for (size_t j = 0; j < D; ++j) {
			r[i] = sub_FxP(r[i], mult_FxP(A[i * D + j], x[j]));
		}
	}

	FxP max_dxj = ZERO;
	do {
		max_dxj = ZERO;
		for (size_t j = 0; j < D; ++j) {
			if (A_norm2[j] == ZERO) continue;
			FxP x_j0 = x[j];
			FxP rho_j = ZERO;
			for (size_t i = 0; i < N; ++i) {
				r[i] = add_FxP(r[i], mult_FxP(A[i * D + j], x[j]));
				rho_j = add_FxP(rho_j, mult_FxP(r[i], A[i * D + j]));
			}
			FxP sign = greater_than_FxP(rho_j, ZERO) ? ONE : MINUS_ONE;
			x[j] = div_FxP(mult_FxP(sign,
																	max_FxP(sub_FxP(abs_FxP(rho_j),
																											lambda),
																						ZERO)),
											A_norm2[j]);
			for (size_t i = 0; i < N; ++i) {
				r[i] = sub_FxP(r[i], mult_FxP(A[i * D + j], x[j]));
			}
			max_dxj = max_FxP(abs_FxP(sub_FxP(x[j], x_j0)), max_dxj);
			max_xj = max_FxP(abs_FxP(x[j]), max_xj);
		}
	} while (greater_than_FxP(div_FxP(max_dxj, max_xj), tol));

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

	FxP *A, *y, *x;
	A = malloc(sizeof(FxP) * N * 2);
	y = malloc(sizeof(FxP) * N);
	x = malloc(sizeof(FxP) * 2);
	x[0] = double_to_FxP(0.5);
	x[1] = double_to_FxP(0.5);
	for(size_t i = 0; i < N; ++i) {
		A[i*2] = double_to_FxP(A_double[i*2]);
		A[i*2+1] = double_to_FxP(A_double[i*2+1]);
		y[i] = double_to_FxP(y_double[i]);
	}
	cd_lasso_FxP(y, A, x, N, 2,
								double_to_FxP(0.1), double_to_FxP(0.0001));
	printf("Estimated parameter (using FxP): (%lf, %lf)\n",
													FxP_to_double(x[0]), FxP_to_double(x[1]));

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
									2.0, 8.0, 32.0, -32.0,
									0.5, 0.25, 0.125, 0.0625,
									1.0e-1, 1.0e-4, 1.0e-8, 1.0e-16,
									1.0e1, 1.0e2, 1.0e4, 1.0e8,
									1.0/3, M_PI, -1.0/3, -1.25,
									-1.0e-1, -1.0e-4, -1.0e-8, -1.0e-16,
									-1.0e1, -1.0e2, -1.0e4, -1.0e8
								 };
	test_FxP_conversion(nums, N);
	const double tol = 1.0e-6;
	test_FxP_arithmetic(nums, N, tol);
	test_FxP_dot_products(20);
	srand(0);
	test_lasso(30);
	return 0;
}
