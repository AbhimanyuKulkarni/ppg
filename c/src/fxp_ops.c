#include "fxp_ops.h"
#include <stdint.h>
#include <stdio.h>
#include <math.h>

FxP int64_to_FxP(int64_t d) {
	return (FxP) ((d << FxP_FRAC_LEN) & (FxP_INT_MASK | FxP_FRAC_MASK));
}

FxP double_to_FxP(double d) {
	union double_bitview u;
	u.dval = d;
	uint64_t double_frac_part = u.uintval & DOUBLE_FRAC_MASK;
	uint64_t double_exp_part = (u.uintval & DOUBLE_EXP_MASK) >> 52;
	uint64_t double_sign_part = (u.uintval & DOUBLE_SIGN_MASK) >> 63;
	
	int64_t fraction;
	if (double_exp_part == 0) {
		// Subnormals.
		// Such values should just become zero though...
		return (FxP) 0;
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
		int frac_shift = (52 - FxP_FRAC_LEN);
		int exp = double_exp_part - 1023;
		int total_shift = frac_shift - exp; // +ve exp -> left shift

		if (exp >= FxP_INT_LEN) {
			// Can't do much here, overflow is guaranteed.
			fprintf(stderr, "Overflow on %lf!\n", d);
		} else if (total_shift > 52) {
			// Underflow. Can't do much...
			fprintf(stderr, "Underflow on %.15lf!\n", d);
			return (FxP) 0;
		}

		if (total_shift > 0) {
			fraction >>= total_shift;
		} else {
			fraction <<= -total_shift;
		}
		if (double_sign_part == 0) {
			return (FxP) fraction;
		} else {
			// perform 2s complement
			// The computer does this for us...
			return (FxP) (-fraction);
		}
	}
}

double FxP_to_double(FxP fxp) {
	union double_bitview u;
	u.uintval = 0;

	int64_t repr = fxp;		// force 64-bit width

	const uint64_t FxP64_SIGN_MASK = ((uint64_t) 1) << 63;
	uint64_t sign_part = repr & FxP64_SIGN_MASK;
	u.uintval |= sign_part;

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
		uint64_t exp = 1023 - FxP_FRAC_LEN + msb_loc;
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

FxP add_FxP(FxP a, FxP b) {
	return a + b;
}

FxP sub_FxP(FxP a, FxP b) {
	return a - b;
}

FxP mult_FxP(FxP a, FxP b) {
	#ifdef PPG_USE_FXP64
		return (FxP) (((((__int128) a) * ((__int128) b)) >> FxP_FRAC_LEN)
										& 0xFFFFFFFFFFFFFFFF);
	#elif defined(PPG_USE_FXP32)
		return (FxP) (((((int64_t) a) * ((int64_t) b)) >> FxP_FRAC_LEN)
										& 0xFFFFFFFF);
	#elif defined(PPG_USE_FXP16)
		return (FxP) (((((int32_t) a) * ((int32_t) b)) >> FxP_FRAC_LEN)
										& 0xFFFF);
	#elif defined(PPG_USE_FXP8)
		return (FxP) (((((int16_t) a) * ((int16_t) b)) >> FxP_FRAC_LEN)
										& 0xFF);
	#else
		return (FxP) 0;	// not implemented
	#endif
}

FxP div_FxP(FxP a, FxP b) {
	#if defined(PPG_USE_FXP64)
		return (FxP) (((((__int128) a) << FxP_FRAC_LEN)
											/ ((__int128) b))
										& 0xFFFFFFFFFFFFFFFF);
	#elif defined(PPG_USE_FXP32)
		return (FxP) (((((int64_t) a) << FxP_FRAC_LEN)
											/ ((int64_t) b))
										& 0xFFFFFFFF);
	#elif defined(PPG_USE_FXP16)
		return (FxP) (((((int32_t) a) << FxP_FRAC_LEN)
											/ ((int32_t) b))
										& 0xFFFF);
	#elif defined(PPG_USE_FXP8)
		return (FxP) (((((int16_t) a) << FxP_FRAC_LEN)
											/ ((int16_t) b))
										& 0xFF);
	#else
		return (FxP) 0; // not implemented
	#endif
}

// signed ints should work the same?
FxP min_FxP(FxP a, FxP b) {
	return a < b ? a : b;
}

FxP max_FxP(FxP a, FxP b) {
	return a > b ? a : b;
}

FxP abs_FxP(FxP x) {
	return x < 0 ? -x : x;
}

bool less_than_FxP(FxP a, FxP b) {
	return sub_FxP(a, b) < 0;
}

bool less_or_eq_FxP(FxP a, FxP b) {
	return sub_FxP(a, b) <= 0;
}

bool greater_than_FxP(FxP a, FxP b) {
	return sub_FxP(a, b) > 0;
}

bool greater_or_eq_FxP(FxP a, FxP b) {
	return sub_FxP(a, b) >= 0;
}

