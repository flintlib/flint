/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_init(nmod_mpoly_t A, const nmod_mpoly_ctx_t ctx)
{
    /* default to MPOLY_MIN_BITS bits per exponent */
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
    A->bits = MPOLY_MIN_BITS;
}

void nmod_mpoly_init3(nmod_mpoly_t A, slong alloc, mp_bitcnt_t bits,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    /* sanitize alloc input */
    alloc = FLINT_MAX(alloc, WORD(0));

    if (alloc != 0)
    {
        A->coeffs = (mp_limb_t *) flint_malloc(alloc*sizeof(mp_limb_t));
        A->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
    } else
    {
        A->coeffs = NULL;
        A->exps = NULL;
    }
    A->alloc = alloc;
    A->length = 0;
    A->bits = bits;
}

void nmod_mpoly_init2(nmod_mpoly_t A, slong alloc, const nmod_mpoly_ctx_t ctx)
{
    /* default to MPOLY_MIN_BITS bits per exponent */
    nmod_mpoly_init3(A, alloc, MPOLY_MIN_BITS, ctx);
}
