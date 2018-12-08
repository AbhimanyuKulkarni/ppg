/*
 * fixed_point_experiments.c
 *
 * Some basic fiddling with fixed point number representations.
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "fixed_point.h"

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

	double sums[32*32];
	double diffs[32*32];
	double prods[32*32];
	double divs[32*32];
	
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
		printf("Testing conversion itself:\n");
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

		printf("----------Testing arithmetic-------------------\n");
		for (size_t j = 0; j < N; ++j) {
			FxP64 fxp_j = double_to_FxP64(nums[j]);

			double sum = FxP64_to_double(add_FxP64(fxp_i, fxp_j));
			double diff = FxP64_to_double(sub_FxP64(fxp_i, fxp_j));
			double prod = FxP64_to_double(mult_FxP64(fxp_i, fxp_j));
			double div = 0.0;
			FxP64 div_actual = 0;
			if (fxp_j != 0) {
				div = FxP64_to_double(div_FxP64(fxp_i, fxp_j));
				div_actual = div_FxP64(fxp_i, fxp_j);
			}
			
			if (sum != sums[i*N+j]) {
				printf("Failed");
			} else {
				printf("OK");
			}
			printf("\t%lf + %lf = %lf, expected %lf\n", nums[i], nums[j],
																									sum, sums[i*N+j]);

			if (diff != diffs[i*N+j]) {
				printf("Failed");
			} else {
				printf("OK");
			}
			printf("\t%lf - %lf = %lf, expected %lf\n", nums[i], nums[j],
																									diff, diffs[i*N+j]);

			if (prod != prods[i*N+j]) {
				printf("Failed");
			} else {
				printf("OK");
			}
			printf("\t%lf * %lf = %lf, expected %lf\n", nums[i], nums[j],
																									prod, prods[i*N+j]);

			if (div != divs[i*N+j]) {
				printf("Failed");
			} else {
				printf("OK");
			}
			printf("\t%lf / %lf = %lf (%016lx), expected %lf\n",
																									nums[i], nums[j],
																									div, div_actual,
																									divs[i*N+j]);
		}
	}

	return 0;
}
