/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef NFLOAT_IMPL_H
#define NFLOAT_IMPL_H

#include "longlong.h"
#include "mpn_extras.h"
#include "nfloat.h"
#include "../mpn_extras/fp_small.h"

/* Whether a result with sign bit sgnbit must be rounded away from zero
   (as opposed to truncated) to respect the directed rounding mode. */
FLINT_FORCE_INLINE int
nfloat_should_round_up(int sgnbit, gr_ctx_t ctx)
{
    if (!NFLOAT_CTX_HAS_DIRECTED_ROUNDING(ctx))
        return 0;
    else
        return sgnbit ^ ((NFLOAT_CTX_FLAGS(ctx) & NFLOAT_RND_CEIL) != 0);
}

/* Checks the exponent of a nonzero finite result. */
#define NFLOAT_CHECK_EXP_RANGE(res, exp, sgnbit, ctx) \
    do { \
        if (FLINT_UNLIKELY((exp) < NFLOAT_MIN_EXP)) \
            return _nfloat_underflow((res), (sgnbit), (ctx)); \
        if (FLINT_UNLIKELY((exp) > NFLOAT_MAX_EXP)) \
            return _nfloat_overflow((res), (sgnbit), (ctx)); \
    } while (0)

/* Increments the normalized n-limb mantissa d of magnitude-truncated
   result, adjusting the exponent on carry out. */
#define NFLOAT_MANT_INCREMENT(d, n, exp) \
    do { \
        if (mpn_add_1((d), (d), (n), 1)) \
        { \
            (d)[(n) - 1] = UWORD(1) << (FLINT_BITS - 1); \
            (exp)++; \
        } \
    } while (0)

/* Copying and zeroing for short operands (mpn_extras.h). */
#define _nfloat_copy_limbs _flint_mpn_copy_small
#define _nfloat_zero_limbs _flint_mpn_zero_small

/* One-limb division kernels (div.c, complex.c): mpn_extras/fp_small.h */
#define NFLOAT_INVERT_LIMB(dinv, d) FLINT_MPN_INVERT_LIMB(dinv, d)
#define _nfloat_div_1_kernel_preinv _flint_mpn_fdiv_1_kernel_preinv
#define _nfloat_div_1_kernel _flint_mpn_fdiv_1_kernel

#endif
