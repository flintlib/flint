/*
    Copyright (C) 2017 William Hart
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_reverse(fmpz_mpoly_t A,
                              const fmpz_mpoly_t B, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    slong Blen = B->length;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    if (A != B)
    {
        fmpz_mpoly_fit_length(A, Blen, ctx);
        fmpz_mpoly_fit_bits(A, B->bits, ctx);
        A->bits = B->bits;
        A->length = B->length;
        for (i = 0; i < Blen; i++)
        {
            fmpz_set(A->coeffs + i,  B->coeffs + Blen - i - 1);
        }
    }
    else
    {
        for (i = 0; i < Blen/2; i++)
        {
            fmpz_swap(A->coeffs + i,  A->coeffs + Blen - i - 1);
        }
    }

    mpoly_reverse(A->exps, B->exps, Blen, N);
}
