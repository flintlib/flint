/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech_mpoly.h"

void fq_zech_mpoly_fit_length_reset_bits(
    fq_zech_mpoly_t A,
    slong len,
    flint_bitcnt_t bits,
    const fq_zech_mpoly_ctx_t ctx)
{
    slong i;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);
    slong new_alloc;

    FLINT_ASSERT(len >= 0);

    if (A->alloc < len)
    {
        new_alloc = FLINT_MAX(len, 2*A->alloc);
        if (A->alloc > 0)
        {
            A->coeffs = (fq_zech_struct *) flint_realloc(A->coeffs, new_alloc*sizeof(fq_zech_struct));
            A->exps = (ulong *) flint_realloc(A->exps, new_alloc*N*sizeof(ulong));
        }
        else
        {
            A->coeffs = (fq_zech_struct *) flint_malloc(new_alloc*sizeof(fq_zech_struct));
            A->exps   = (ulong *) flint_malloc(new_alloc*N*sizeof(ulong));
        }
        for (i = A->alloc; i < new_alloc; i++)
            fq_zech_init(A->coeffs + i, ctx->fqctx);
        A->alloc = new_alloc;
    }
    else if (A->bits < bits)
    {
        if (A->alloc > 0)
            A->exps = (ulong *) flint_realloc(A->exps, A->alloc*N*sizeof(ulong));
    }

    A->bits = bits;
}

