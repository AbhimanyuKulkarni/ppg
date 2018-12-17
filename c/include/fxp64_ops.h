#ifndef PPG_FXP64_OPS_H
#define PPG_FXP64_OPS_H

#include "fixed_point.h"
#include <math.h>

// convert x to/from right type
#define FLOAT64_TO_NUM(x) (double_to_FxP64(x))
#define UINT_TO_NUM(x) ((((int64_t) x) << FxP64_FRAC_LEN) & 0xFFFFFFFFFFFFFFFF)
#define NUM_TO_FLOAT64(x) (FxP64_to_double(x))

// arithmetic and logical
#define ADD(a, b) add_FxP64(a, b)
#define SUB(a, b) sub_FxP64(a, b)
#define MULT(a, b) mult_FxP64(a, b)
#define DIV(a, b) div_FxP64(a, b)
#define GREATER_THAN(a, b) (a > b)
#define GREATER_OR_EQ(a, b) (a >= b)
#define LESS_THAN(a, b) (a < b)
#define LESS_OR_EQ(a, b) (a <= b)
#define EQUAL(a, b) (a == b)

// more math
// see http://graphics.stanford.edu/~seander/bithacks.html#IntegerAbs
#define FABS(x) (LESS_THAN(x, 0) ? -x : x)
#define SQRT(x) (double_to_FxP64(sqrt(FxP64_to_double(x))))
#define COS(x) (double_to_FxP64(cos(FxP64_to_double(x))))
#endif
