/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_init(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);

    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fq_nmod_mpoly_init3(fq_nmod_mpoly_t A, slong alloc, flint_bitcnt_t bits,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    if (alloc != 0)
    {
        slong i;
        A->coeffs = (fq_nmod_struct *) flint_malloc(alloc*sizeof(fq_nmod_struct));
        A->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
        for (i = 0; i < alloc; i++)
            fq_nmod_init(A->coeffs + i, ctx->fqctx);
    }
    else
    {
        A->coeffs = NULL;
        A->exps = NULL;
    }
    A->alloc = alloc;
    A->length = 0;
    A->bits = bits;
}

void fq_nmod_mpoly_init2(fq_nmod_mpoly_t A, slong alloc,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);
    fq_nmod_mpoly_init3(A, alloc, bits, ctx);
}
