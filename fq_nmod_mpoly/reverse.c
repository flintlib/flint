/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_reverse(fq_nmod_mpoly_t A,
                        const fq_nmod_mpoly_t B, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong Blen = B->length;
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    if (A != B)
    {
        fq_nmod_mpoly_fit_length(A, Blen, ctx);
        fq_nmod_mpoly_fit_bits(A, B->bits, ctx);
        A->bits = B->bits;
        A->length = B->length;
        for (i = 0; i < Blen; i++)
        {
            fq_nmod_set(A->coeffs + i,  B->coeffs + Blen - i - 1, ctx->fqctx);
        }
    }
    else
    {
        for (i = 0; i < Blen/2; i++)
        {
            fq_nmod_swap(A->coeffs + i,  A->coeffs + Blen - i - 1, ctx->fqctx);
        }
    }

    mpoly_reverse(A->exps, B->exps, Blen, N);
}
