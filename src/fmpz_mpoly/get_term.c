/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly.h"

void fmpz_mpoly_get_term(fmpz_mpoly_t M, const fmpz_mpoly_t A,
                                           slong i, const fmpz_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t bits = A->bits;

    if ((ulong) i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpz_mpoly_get_term");
    }

    fmpz_mpoly_fit_length(M, WORD(1), ctx);
    fmpz_mpoly_fit_bits(M, bits, ctx);
    M->bits = bits;

    N = mpoly_words_per_exp(bits, ctx->minfo);

    mpoly_monomial_set(M->exps + N*0, A->exps + N*i, N);
    fmpz_set(M->coeffs + 0, A->coeffs + i);
    _fmpz_mpoly_set_length(M, 1, ctx);
}
