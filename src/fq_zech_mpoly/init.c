/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void fq_zech_mpoly_init(fq_zech_mpoly_t A, const fq_zech_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);

    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = bits;
}

void fq_zech_mpoly_init3(fq_zech_mpoly_t A, slong alloc, flint_bitcnt_t bits,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    FLINT_ASSERT((MPOLY_MIN_BITS <= bits && bits <= FLINT_BITS) ||
                 (bits % FLINT_BITS) == 0);

    if (alloc > 0)
    {
        slong i;
        A->coeffs = (fq_zech_struct *) flint_malloc(alloc*sizeof(fq_zech_struct));
        A->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
        for (i = 0; i < alloc; i++)
            fq_zech_init(A->coeffs + i, ctx->fqctx);
        A->alloc = alloc;
    }
    else
    {
        A->coeffs = NULL;
        A->exps = NULL;
        A->alloc = 0;
    }
    A->length = 0;
    A->bits = bits;
}

void fq_zech_mpoly_init2(fq_zech_mpoly_t A, slong alloc,
                                                 const fq_zech_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);
    fq_zech_mpoly_init3(A, alloc, bits, ctx);
}
