/*
    Copyright (C) 2016 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_monomial_fmpz(fmpz ** exps, const fmpz_mpoly_t poly, 
                                           slong n, const fmpz_mpoly_ctx_t ctx)
{
    slong i, N;
    fmpz * temp;
    TMP_INIT;

    TMP_START;
    temp = (fmpz *) TMP_ALLOC(ctx->minfo->nvars*sizeof(fmpz));
    for (i = 0; i < ctx->minfo->nvars; i++)
        fmpz_init(temp + i);

    N = mpoly_words_per_exp(poly->bits, ctx->minfo);
    mpoly_get_monomial_fmpz(temp, poly->exps + N*n, poly->bits, ctx->minfo);

    for (i = 0; i < ctx->minfo->nvars; i++)
    {
        fmpz_swap(temp + i, exps[i]);
        fmpz_clear(temp + i);
    }
    TMP_END;
}
