/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_init(nmod_mpoly_t poly, const nmod_mpoly_ctx_t ctx)
{
    /* default to at least MPOLY_MIN_BITS bits per exponent */
    slong bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);

    poly->coeffs = NULL;
    poly->exps = NULL;
    poly->alloc = 0;
    poly->length = 0;
    poly->bits = bits;
}

void nmod_mpoly_init3(nmod_mpoly_t poly, slong alloc, mp_bitcnt_t bits,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    if (alloc != 0)
    {
        poly->coeffs = (mp_limb_t *) flint_malloc(alloc*sizeof(mp_limb_t));
        poly->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
    } else
    {
        poly->coeffs = NULL;
        poly->exps = NULL;
    }
    poly->alloc = alloc;
    poly->length = 0;
    poly->bits = bits;
}

void nmod_mpoly_init2(nmod_mpoly_t poly, slong alloc, const nmod_mpoly_ctx_t ctx)
{
    /* default to at least MPOLY_MIN_BITS bits per exponent */
    slong bits = mpoly_fix_bits(MPOLY_MIN_BITS, ctx->minfo);
    nmod_mpoly_init3(poly, alloc, bits, ctx);
}

