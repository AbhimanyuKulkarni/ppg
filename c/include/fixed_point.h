/*
 * fixed_point.h
 *
 * See https://en.wikipedia.org/wiki/Double-precision_floating-point_format
 * for where I have got these definitions.
 *
 */
#ifndef PPG_FxP_H
#define PPG_FxP_H

#include <stdint.h>

typedef int64_t FxP64;
typedef int32_t FxP32;
typedef int16_t FxP16;
typedef int8_t FxP8;

#define DOUBLE_SIGN_MASK ((int64_t) (((int64_t) 1) << 63))
#define DOUBLE_EXP_MASK ((int64_t) (((int64_t) 0x7FF) << 52))
#define DOUBLE_FRAC_MASK ((int64_t) ((((int64_t) 1) << 52) - 1))
#define DOUBLE_ONE_FLAG ((int64_t) 1 << 52)

#define FxP64_INT_LEN 16
#define FxP64_FRAC_LEN 48
#define FxP64_INT_MASK ((int64_t) ((((int64_t) 1) << FxP64_INT_LEN) - 1) << FxP64_FRAC_LEN)
#define FxP64_FRAC_MASK ((int64_t) ((((int64_t) 1) << FxP64_FRAC_LEN) - 1))
#define FxP64_SIGN_MASK ((int64_t) (((int64_t) 1) << 63))

union double_bitview {
	double dval;
	uint64_t uintval;
};

FxP64 int64_to_FxP64(int64_t d);
FxP64 double_to_FxP64(double d);
double FxP64_to_double(FxP64 fxp);

// Will be used by fxp64_ops.h
FxP64 add_FxP64(FxP64 a, FxP64 b);
FxP64 sub_FxP64(FxP64 a, FxP64 b);
FxP64 mult_FxP64(FxP64 a, FxP64 b);
FxP64 div_FxP64(FxP64 a, FxP64 b);
#endif
