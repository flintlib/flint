/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* this file DOES NOT need to change with new orderings */

/*
    Get a user exponent "exp"' such that it can be packed into "exp_bits" bits.
    The count "exp_bits" includes the extra bits for the sign.
*/
void mpoly_monomial_randbits_fmpz(fmpz * exp, flint_rand_t state,
                                  flint_bitcnt_t exp_bits, const mpoly_ctx_t mctx)
{
    slong j;
    flint_bitcnt_t newbits = exp_bits;

    while (newbits != (flint_bitcnt_t)(0))
    {
        for (j = 0; j < mctx->nvars; j++) {
            fmpz_randtest_unsigned(exp + j, state, newbits);
        }

        if (mpoly_exp_bits_required_ffmpz(exp, mctx) <= exp_bits)
            return;

        newbits--;
    }

    for (j = 0; j < mctx->nvars; j++) {
        fmpz_zero(exp + j);
    }
}
