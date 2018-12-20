#ifndef PPG_FXP_OPS_H
#define PPG_FXP_OPS_H
/*
 * fxp_ops.h
 *
 * See https://en.wikipedia.org/wiki/Double-precision_floating-point_format
 * for where I have got these definitions.
 */

#include <math.h>
#include <stdint.h>
#include <stdbool.h>

// The file which includes this should define PPG_FRAC already.
// Right now that file is setup.h
#ifdef PPG_USE_FXP8
	typedef int8_t FxP;
	#define FxP_INT_LEN 4
	#define FxP_FRAC_LEN 4
	#define FxP_INT_MASK ((int8_t) ((((int8_t) 1) << FxP_INT_LEN) - 1) << FxP_FRAC_LEN)
	#define FxP_FRAC_MASK ((int8_t) ((((int8_t) 1) << FxP_FRAC_LEN) - 1))
	#define FxP_SIGN_MASK ((int8_t) (((int8_t) 1) << 7))
#elif defined(PPG_USE_FXP16)
	typedef int16_t FxP;
	#define FxP_INT_LEN 8
	#define FxP_FRAC_LEN 8
	#define FxP_INT_MASK ((int16_t) ((((int16_t) 1) << FxP_INT_LEN) - 1) << FxP_FRAC_LEN)
	#define FxP_FRAC_MASK ((int16_t) ((((int16_t) 1) << FxP_FRAC_LEN) - 1))
	#define FxP_SIGN_MASK ((int16_t) (((int16_t) 1) << 15))
#elif defined(PPG_USE_FXP32)
	typedef int32_t FxP;
	#define FxP_INT_LEN 16
	#define FxP_FRAC_LEN 16
	#define FxP_INT_MASK ((int32_t) ((((int32_t) 1) << FxP_INT_LEN) - 1) << FxP_FRAC_LEN)
	#define FxP_FRAC_MASK ((int32_t) ((((int32_t) 1) << FxP_FRAC_LEN) - 1))
	#define FxP_SIGN_MASK ((int32_t) (((int32_t) 1) << 31))
#elif defined(PPG_USE_FXP64)
	typedef int64_t FxP;
	#define FxP_INT_LEN 32
	#define FxP_FRAC_LEN 32
	#define FxP_INT_MASK ((int64_t) ((((int64_t) 1) << FxP_INT_LEN) - 1) << FxP_FRAC_LEN)
	#define FxP_FRAC_MASK ((int64_t) ((((int64_t) 1) << FxP_FRAC_LEN) - 1))
	#define FxP_SIGN_MASK ((int64_t) (((int64_t) 1) << 63))
#else	// default
	typedef int64_t FxP;
	#define FxP_INT_LEN 32
	#define FxP_FRAC_LEN 32
	#define FxP_INT_MASK ((int64_t) ((((int64_t) 1) << FxP_INT_LEN) - 1) << FxP_FRAC_LEN)
	#define FxP_FRAC_MASK ((int64_t) ((((int64_t) 1) << FxP_FRAC_LEN) - 1))
	#define FxP_SIGN_MASK ((int64_t) (((int64_t) 1) << 63))
#endif

#define DOUBLE_SIGN_MASK ((int64_t) (((int64_t) 1) << 63))
#define DOUBLE_EXP_MASK ((int64_t) (((int64_t) 0x7FF) << 52))
#define DOUBLE_FRAC_MASK ((int64_t) ((((int64_t) 1) << 52) - 1))
#define DOUBLE_ONE_FLAG ((int64_t) 1 << 52)

union double_bitview {
	double dval;
	uint64_t uintval;
};

// convert x to/from right type
#define FLOAT64_TO_FRAC(x) (double_to_FxP(x))
#define INT64_TO_FRAC(x) (int64_to_FxP(x))
#define FRAC_TO_FLOAT64(x) (FxP_to_double(x))

// arithmetic and logical
#define ADD(a, b) add_FxP(a, b)
#define SUB(a, b) sub_FxP(a, b)
#define MULT(a, b) mult_FxP(a, b)
#define DIV(a, b) div_FxP(a, b)
#define GREATER_THAN(a, b) (greater_than_FxP(a, b))
#define GREATER_OR_EQ(a, b) (greater_or_eq_FxP(a, b))
#define LESS_THAN(a, b) (less_than_FxP(a, b))
#define LESS_OR_EQ(a, b) (less_or_eq_FxP(a, b))
#define EQUAL(a, b) (a == b)	// just checking bits directly pretty much

#define MAX(a, b) (max_FxP(a, b))
#define MIN(a, b) (min_FxP(a, b))

// more math
#define FABS(x) (abs_FxP(x))

FxP int64_to_FxP(int64_t d);
FxP double_to_FxP(double d);
double FxP_to_double(FxP fxp);

FxP add_FxP(FxP a, FxP b);
FxP sub_FxP(FxP a, FxP b);
FxP mult_FxP(FxP a, FxP b);
FxP div_FxP(FxP a, FxP b);
FxP min_FxP(FxP a, FxP b);
FxP max_FxP(FxP a, FxP b);
FxP abs_FxP(FxP x);
bool less_than_FxP(FxP a, FxP b);
bool less_or_eq_FxP(FxP a, FxP b);
bool greater_than_FxP(FxP a, FxP b);
bool greater_or_eq_FxP(FxP a, FxP b);
#endif
