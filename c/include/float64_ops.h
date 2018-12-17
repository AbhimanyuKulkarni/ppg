#ifndef PPG_FLOAT64_OPS_H
#define PPG_FLOAT64_OPS_H

#include <math.h>

// convert (double) x to/from right type
// double is already the right type if we're using double
#define FLOAT64_TO_NUM(x) (x)
#define UINT_TO_NUM(x) ((double) x)
#define NUM_TO_FLOAT64(x) (x)

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

// more math
#define FABS(x) fabs(x)
#define SQRT(x) sqrt(x)
#define COS(x) cos(x)
#endif