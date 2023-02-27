/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"


void mpoly_monomials_deflation(fmpz * shift, fmpz * stride,
                        const ulong * Aexps, flint_bitcnt_t Abits, slong Alength,
                                                        const mpoly_ctx_t mctx)
{
    slong i, j;
    slong NA;
    slong nvars = mctx->nvars;
    fmpz_t d;
    fmpz * exps;
    TMP_INIT;

    for (j = 0; j < nvars; j++)
        fmpz_zero(stride + j);

    if (Alength == 0)
    {
        /* undefined case */
        for (j = 0; j < nvars; j++)
            fmpz_zero(shift + j);
        return;
    }

    TMP_START;

    exps = (fmpz *) TMP_ALLOC(nvars*sizeof(fmpz));
    for (j = 0; j < nvars; j++)
        fmpz_init(exps + j);
    
    fmpz_init(d);

    NA = mpoly_words_per_exp(Abits, mctx);

    i = Alength - 1;
    mpoly_get_monomial_ffmpz(shift, Aexps + NA*i, Abits, mctx);

    for (i--; i >= 0; i--)
    {
        mpoly_get_monomial_ffmpz(exps, Aexps + NA*i, Abits, mctx);
        for (j = 0; j < nvars; j++)
        {
            fmpz_sub(d, exps + j, shift + j);
            fmpz_gcd(stride + j, stride + j, d);
            if (fmpz_sgn(d) < 0)
                fmpz_swap(shift + j, exps + j);
        }
    }

    for (j = 0; j < nvars; j++)
        fmpz_clear(exps + j);

    fmpz_clear(d);

    TMP_END;
}

