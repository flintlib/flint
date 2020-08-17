/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_get_term(fmpq_mpoly_t M, const fmpq_mpoly_t A,
                                           slong i, const fmpq_mpoly_ctx_t ctx)
{
    slong N;
    flint_bitcnt_t bits = A->zpoly->bits;

    if ((ulong) i >= (ulong) A->zpoly->length)
    {
        flint_throw(FLINT_ERROR, "Index out of range in fmpq_mpoly_get_term");
    }

    fmpz_mpoly_fit_length(M->zpoly, WORD(1), ctx->zctx);
    fmpz_mpoly_fit_bits(M->zpoly, bits, ctx->zctx);
    M->zpoly->bits = bits;

    N = mpoly_words_per_exp(bits, ctx->zctx->minfo);

    mpoly_monomial_set(M->zpoly->exps + N*0, A->zpoly->exps + N*i, N);
    fmpq_mul_fmpz(M->content, A->content, A->zpoly->coeffs + i);
    fmpz_one(M->zpoly->coeffs + 0);
    _fmpz_mpoly_set_length(M->zpoly, 1, ctx->zctx);
}
