/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_reverse(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);
    slong i;
    slong Blen = B->length;

    if (A != B)
    {
        fq_nmod_mpoly_fit_length_reset_bits(A, Blen, B->bits, ctx);
        A->length = B->length;
        for (i = 0; i < Blen; i++)
            _n_fq_set(A->coeffs + d*i, B->coeffs + d*(Blen - i - 1), d);
    }
    else
    {
        for (i = 0; i < Blen/2; i++)
            _n_fq_swap(A->coeffs + d*i,  A->coeffs + d*(Blen - i - 1), d);
    }

    mpoly_reverse(A->exps, B->exps, Blen, N);
}
