/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/* !!! this file DOES need to change with new orderings */

/*
    compute number of bits required to store user_exp in packed format
    the returned number of bits includes space for a zero'd signed bit
*/
flint_bitcnt_t mpoly_exp_bits_required_ui(const ulong * user_exp,
                                                        const mpoly_ctx_t mctx)
{
    slong i, nfields = mctx->nfields;
    ulong max = 0;

    if (mctx->deg)
    {
        for (i = 0; i < nfields - 1; i++)
        {
            max += user_exp[i];
            if (max < user_exp[i])
                return 2*FLINT_BITS;
        }
    }
    else
    {
        for (i = 0; i < nfields; i++)
        {
            max |= user_exp[i];
        }
    }

    return 1 + FLINT_BIT_COUNT(max);
}

/*
    compute number of bits required to store user_exp in packed format
    the returned number of bits includes space for a zero'd signed bit
*/
flint_bitcnt_t mpoly_exp_bits_required_ffmpz(const fmpz * user_exp,
                                                        const mpoly_ctx_t mctx)
{
    slong i, nvars = mctx->nvars;
    flint_bitcnt_t exp_bits;

    if (mctx->deg)
    {
        fmpz_t deg;
        fmpz_init(deg);
        for (i = 0; i < nvars; i++)
        {
            fmpz_add(deg, deg, user_exp + i);
        }
        exp_bits = 1 + fmpz_bits(deg);
        fmpz_clear(deg);
    }
    else
    {
        exp_bits = 0;
        for (i = 0; i < nvars; i++)
        {
            flint_bitcnt_t this_bits = fmpz_bits(user_exp + i);
            exp_bits = FLINT_MAX(exp_bits, this_bits);
        }
        exp_bits += 1;
    }

    return exp_bits;
}

/*
    compute number of bits required to store user_exp in packed format
    the returned number of bits includes space for a zero'd signed bit
*/
flint_bitcnt_t mpoly_exp_bits_required_pfmpz(fmpz * const * user_exp,
                                                        const mpoly_ctx_t mctx)
{
    slong i, nvars = mctx->nvars;
    flint_bitcnt_t exp_bits;

    if (mctx->deg)
    {
        fmpz_t deg;
        fmpz_init(deg);
        for (i = 0; i < nvars; i++)
        {
            fmpz_add(deg, deg, user_exp[i]);
        }
        exp_bits = 1 + fmpz_bits(deg);
        fmpz_clear(deg);
    }
    else
    {
        exp_bits = 0;
        for (i = 0; i < nvars; i++)
        {
            flint_bitcnt_t this_bits = fmpz_bits(user_exp[i]);
            exp_bits = FLINT_MAX(exp_bits, this_bits);
        }
        exp_bits += 1;
    }

    return exp_bits;
}
