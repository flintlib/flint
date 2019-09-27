/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_deflate(fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
        const fmpz * shift, const fmpz * stride, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i;
    slong Abits;

    if (B->length == 0)
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    /* quick and safe bound on bits required */
    Abits = B->bits;

    if (A == B)
    {
        slong NA = mpoly_words_per_exp(Abits, ctx->minfo);
        ulong * texps = flint_malloc(NA*A->alloc*sizeof(ulong));
        mpoly_monomials_deflate(texps, Abits, B->exps, B->bits, B->length,
                                                    shift, stride, ctx->minfo);
        flint_free(A->exps);
        A->exps = texps;
        A->bits = Abits;
    }
    else
    {
        fq_nmod_mpoly_fit_length(A, B->length, ctx);
        fq_nmod_mpoly_fit_bits(A, Abits, ctx);
        A->bits = Abits;
        for (i = 0; i < B->length; i++)
            fq_nmod_set(A->coeffs + i, B->coeffs + i, ctx->fqctx);
        mpoly_monomials_deflate(A->exps, Abits, B->exps, B->bits, B->length,
                                                    shift, stride, ctx->minfo);
        _fq_nmod_mpoly_set_length(A, B->length, ctx);
    }

    if (ctx->minfo->ord != ORD_LEX)
    {
        fq_nmod_mpoly_sort_terms(A, ctx);
    }
}
