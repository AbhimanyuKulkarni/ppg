#include "fixed_point.h"
#include <stdint.h>
#include <stdio.h>
#include <math.h>

FxP64 int64_to_FxP64(int64_t d) {
	return d << FxP64_FRAC_LEN;
}

FxP64 double_to_FxP64(double d) {
	union double_bitview u;
	u.dval = d;
	uint64_t double_frac_part = u.uintval & DOUBLE_FRAC_MASK;
	uint64_t double_exp_part = (u.uintval & DOUBLE_EXP_MASK) >> 52;
	uint64_t double_sign_part = (u.uintval & DOUBLE_SIGN_MASK) >> 63;
	
	int64_t fraction;
	if (double_exp_part == 0) {
		// Subnormals.
		// Such values should just become zero though...
		return (FxP64) 0;
	} else {
		// Number is -1^sign * 2^(e-1023) * 1.<fraction>
		fraction = double_frac_part | DOUBLE_ONE_FLAG;

		// => left/right shift 1.<fraction> by (e - 1023) bits, and take the
		// MSBs of the result?
		// The fraction is 53 bits right now. These 53 bits have to be
		// brought right first (if FxP64_FRAC_LEN < 52) to fit in the
		// fractional part (so right shift by (52 - FxP64_FRAC_LEN)).
		// But then we left shift by double_exp_part, so let's combine the
		// two shifts and do them together.
		// This is not a big deal with 64 bits, but should make a
		// substantial difference when we extend this to shorter
		// representations (32-, 16-, and especially 8-bit).
		int frac_shift = (52 - FxP64_FRAC_LEN);
		int exp = double_exp_part - 1023;
		int total_shift = frac_shift - exp; // +ve exp -> left shift

		if (exp >= FxP64_INT_LEN) {
			// Can't do much here, overflow is guaranteed.
			fprintf(stderr, "Overflow on %lf!\n", d);
		} else if (total_shift > 52) {
			// Underflow. Can't do much...
			fprintf(stderr, "Underflow on %.15lf!\n", d);
		}

		if (total_shift > 0) {
			fraction >>= total_shift;
		} else {
			fraction <<= -total_shift;
		}
		if (double_sign_part == 0) {
			return fraction;
		} else {
			// perform 2s complement
			// The computer does this for us...
			return -fraction;
		}
	}
}

double FxP64_to_double(FxP64 fxp) {
	union double_bitview u;
	u.uintval = 0;

	uint64_t sign_part = fxp & FxP64_SIGN_MASK;
	u.uintval |= sign_part;
	
	uint64_t repr = fxp;
	if (sign_part) {
		// re-do 2s complement to get magnitude.
		// But the computer can do this for us automatically...
		repr = -repr;
	}

	if (repr != 0) {
		// Find MSB
		uint8_t msb_loc = (uint8_t) floor(log2(repr));
		// Exponent is related to MSB:
		// If msb = 0, exponent is -FxP64_FRAC_LEN
		// if msb = FxP64_FRAC_LEN, exponent is 0
		// if msb = 62, exponent is FxP64_INT_LEN-2
		//							already max was FxP64_INT_LEN-1
		//							if unsigned, but with the last bit being sign
		//							we should have max exponent being FxP64_INT_LEN-2
		//							e.g. for 16-48
		//								msb 0 => 2^-48 is the max term
		//										  => exponent is -48 then
		//								msb 48 => number is already 1.<fraction>
		//								so if msb 62
		//									we should have max term 2^14
		uint64_t exp = 1023 - FxP64_FRAC_LEN + msb_loc;
		exp <<= 52;
		u.uintval |= exp;

		// Now for the magnitude
		// Number is 1.<something> from MSB onwards.
		// Which is what we want...
		// We need to get the MSB to location 52, and then take locations
		// 0-51 for the fraction part in the number
		if (msb_loc < 52) {
			repr <<= (52 - msb_loc);
		} else {
			repr >>= (msb_loc - 52);
		}
		
		u.uintval |= (repr & DOUBLE_FRAC_MASK);
	}
	return u.dval;
}

FxP64 add_FxP64(FxP64 a, FxP64 b) {
	return a + b;
}

FxP64 sub_FxP64(FxP64 a, FxP64 b) {
	return a - b;
}

FxP64 mult_FxP64(FxP64 a, FxP64 b) {
	return (FxP64) (((((__int128) a) * ((__int128) b)) >> FxP64_FRAC_LEN)
									& 0xFFFFFFFFFFFFFFFF);
}

FxP64 div_FxP64(FxP64 a, FxP64 b) {
	return (FxP64) (((((__int128) a) << FxP64_FRAC_LEN) 
										/ ((__int128) b))
									& 0xFFFFFFFFFFFFFFFF);
}
