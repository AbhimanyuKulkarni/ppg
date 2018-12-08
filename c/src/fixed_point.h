#ifndef PPG_FxP_H
#define PPG_FxP_H

#include <stdint.h>

typedef uint64_t FxP64;
typedef uint32_t FxP32;
typedef uint16_t FxP16;
typedef uint8_t FxP8;

#define DOUBLE_SIGN_MASK ((uint64_t) (((uint64_t) 1) << 63))
#define DOUBLE_EXP_MASK ((uint64_t) (((uint64_t) 0x7FF) << 52))
#define DOUBLE_FRAC_MASK ((uint64_t) ((((uint64_t) 1) << 52) - 1))
#define DOUBLE_ONE_FLAG ((uint64_t) 1 << 52)

#define FxP64_INT_LEN 16
#define FxP64_FRAC_LEN 48
#define FxP64_INT_MASK ((uint64_t) ((((uint64_t) 1) << FxP64_INT_LEN) - 1) << FxP64_FRAC_LEN)
#define FxP64_FRAC_MASK ((uint64_t) ((((uint64_t) 1) << FxP64_FRAC_LEN) - 1))
#define FxP64_SIGN_MASK ((uint64_t) (((uint64_t) 1) << 63))

union double_bitview {
	double dval;
	uint64_t uintval;
};

FxP64 double_to_FxP64(double d);
double FxP64_to_double(FxP64 fxp);
#endif
