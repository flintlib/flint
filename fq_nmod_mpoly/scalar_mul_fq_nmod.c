/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_scalar_mul_fq_nmod(fq_nmod_mpoly_t A,
     const fq_nmod_mpoly_t B, const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;

    if (fq_nmod_is_zero(c, ctx->fqctx))
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    if (A == B)
    {
        if (fq_nmod_is_one(c, ctx->fqctx))
            return;
    }
    else
    {
        slong N;

        fq_nmod_mpoly_fit_length(A, B->length, ctx);
        fq_nmod_mpoly_fit_bits(A, B->bits, ctx);
        A->length = B->length;
        A->bits = B->bits;

        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        mpoly_copy_monomials(A->exps, B->exps, B->length, N);
        if (fq_nmod_is_one(c, ctx->fqctx))
        {
            for (i = 0; i < B->length; i++)
            {
                fq_nmod_set(A->coeffs + i, B->coeffs + i, ctx->fqctx);
            }
            return;
        }
    }

    for (i = 0; i < B->length; i++)
    {
        fq_nmod_mul(A->coeffs + i, B->coeffs + i, c, ctx->fqctx);
    }
}
