/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_init3(
    fq_nmod_mpoly_t A,
    slong alloc,
    flint_bitcnt_t bits,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    if (alloc > 0)
    {
        A->coeffs_alloc = d*alloc;
        A->coeffs = FLINT_ARRAY_ALLOC(A->coeffs_alloc, mp_limb_t);
        A->exps_alloc = N*alloc;
        A->exps = FLINT_ARRAY_ALLOC(A->exps_alloc, ulong);
    }
    else
    {
        A->coeffs = NULL;
        A->exps = NULL;
        A->coeffs_alloc = 0;
        A->exps_alloc = 0;
    }
    A->length = 0;
    A->bits = bits;
}

void fq_nmod_mpoly_init2(
    fq_nmod_mpoly_t A,
    slong alloc,
    const fq_nmod_mpoly_ctx_t ctx)
{
    flint_bitcnt_t bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);
    fq_nmod_mpoly_init3(A, alloc, bits, ctx);
}

