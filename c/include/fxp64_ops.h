#ifndef PPG_FXP64_OPS_H
#define PPG_FXP64_OPS_H

#include "fixed_point.h"
#include <math.h>

// convert x to/from right type
#define FLOAT64_TO_FRAC(x) (double_to_FxP64(x))
#define INT64_TO_FRAC(x) (int64_to_FxP64(x))
#define FRAC_TO_FLOAT64(x) (FxP64_to_double(x))

// arithmetic and logical
#define ADD(a, b) add_FxP64(a, b)
#define SUB(a, b) sub_FxP64(a, b)
#define MULT(a, b) mult_FxP64(a, b)
#define DIV(a, b) div_FxP64(a, b)
#define GREATER_THAN(a, b) (greater_than_FxP64(a, b))
#define GREATER_OR_EQ(a, b) (greater_or_eq_FxP64(a, b))
#define LESS_THAN(a, b) (less_than_FxP64(a, b))
#define LESS_OR_EQ(a, b) (less_or_eq_FxP64(a, b))
#define EQUAL(a, b) (a == b)	// just checking bits directly pretty much

#define MAX(a, b) (max_FxP64(a, b))

// more math
#define FABS(x) (abs_FxP64(x))
#endif
