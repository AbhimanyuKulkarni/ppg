/*
 * fixed_point_experiments.c
 *
 * Some basic fiddling with fixed point number representations.
 */

#include <stdio.h>
#include <stdint.h>
#include <math.h>

typedef int32_t FxP;

#define FxP_32BIT_MASK 0xFFFFFFFF
#define FxP_INT_MASK 0xFFFFFF00
#define FxP_FRAC_MASK 0x000000FF

#define ADD_FXP(a, b) (a + b)
#define SUB_FXP(a, b) (a - b)
#define MULTIPLY_FXP(a, b) ((FxP) (((int64_t) a * (int64_t) b) & FxP_32BIT_MASK))
#define DIVIDE_FXP(a, b) ((FxP) ((((int64_t) a << 32) / (int64_t) b) >> 32))

FxP double_to_FxP(double f) {
	double intpart_f, fracpart_f;
	// fracpart_f = fabs(modf(f, &intpart_f));
	intpart_f = floor(f);
	fracpart_f = f - intpart_f;
	int32_t intpart = (int32_t) intpart_f;
	int32_t fracpart = (int32_t) (fracpart_f * (1 << 8));
	return (FxP_INT_MASK & (intpart << 8)) + (FxP_FRAC_MASK & fracpart);
}

int main(int argc, char *argv[]) {
	FxP fxp1 = double_to_FxP(1.5);
	FxP fxp2 = double_to_FxP(1.5);
	FxP fxp3 = double_to_FxP(3.375);
	printf("%d * %d = %d, should be %d\n", fxp1, fxp2, MULTIPLY_FXP(fxp1, fxp2), fxp3);
	return 0;
}
