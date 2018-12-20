#ifndef PPG_FLOAT_OPS_H
#define PPG_FLOAT_OPS_H

#include <math.h>

// The file which includes this should define PPG_FRAC already.
// Right now that file is setup.h

// convert (double) x to/from right type
// With floating point there's not much to do
// 2 cases: float or double; both cases just need casts
// Need to be careful with ranges of values.
#define FLOAT64_TO_FRAC(x) ((PPG_FRAC) x)
#define INT64_TO_FRAC(x) ((PPG_FRAC) x)
#define FRAC_TO_FLOAT64(x) ((double) x)

// arithmetic and logical
#define ADD(a, b) (a + b)
#define SUB(a, b) (a - b)
#define MULT(a, b) (a * b)
#define DIV(a, b) (a / b)
#define GREATER_THAN(a, b) (a > b)
#define GREATER_OR_EQ(a, b) (a >= b)
#define LESS_THAN(a, b) (a < b)
#define LESS_OR_EQ(a, b) (a <= b)
#define EQUAL(a, b) (a == b)

#define MAX(a, b) (a > b ? a : b)
#define MIN(a, b) (a < b ? a : b)

// more math
#define FABS(x) fabs(x)
#endif
