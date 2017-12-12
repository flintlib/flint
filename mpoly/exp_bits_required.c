/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "mpoly.h"

/*
    compute number of bits required to store user_exp in packed format
    the returned number of bits includes space for a zero'd signed bit
    a return value of > FLINT_BITS indicates an error (it doesn't fit)
*/
slong mpoly_exp_bits_required_ui(const ulong * user_exp, const mpoly_ctx_t mctx)
{
    int deg = mctx->deg;
    slong i, exp_bits, nfields = mctx->nfields;
    ulong max = 0;

    if (deg)
    {
        for (i = 0; i < nfields - 1; i++)
        {
            max += user_exp[i];
            if (max < user_exp[i])
                return 2*FLINT_BITS;
        }
    } else
    {
        for (i = 0; i < nfields; i++)
        {
            if (max < user_exp[i])
                max = user_exp[i];
        }
    }

    exp_bits = FLINT_MAX(WORD(8), FLINT_BIT_COUNT(max) + 1);
    return exp_bits;
}

/*
    compute number of bits required to store user_exp in packed format
    the returned number of bits includes space for a zero'd signed bit
*/
mp_bitcnt_t mpoly_exp_bits_required_fmpz(const fmpz * user_exp,
                                                        const mpoly_ctx_t mctx)
{
    int deg = mctx->deg;
    slong i, exp_bits, nvars = mctx->nvars;
    fmpz_t max;
    fmpz_init(max);
    fmpz_zero(max);

    if (deg)
    {
        for (i = 0; i < nvars; i++)
            fmpz_add(max, max, user_exp + i);

    } else
    {
        for (i = 0; i < nvars; i++)
            if (fmpz_cmp(max, user_exp + i) < 0)
                fmpz_set(max, user_exp + i);
    }

    exp_bits = fmpz_bits(max) + 1;

    fmpz_clear(max);

    return exp_bits;
}
