/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_neg(nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    slong i, N;

    if (A != B)
    {
        N = mpoly_words_per_exp(B->bits, ctx->minfo);
        nmod_mpoly_fit_length(A, B->length, ctx);
        nmod_mpoly_fit_bits(A, B->bits, ctx);
        A->bits = B->bits;
        mpn_copyi(A->exps, B->exps, N*B->length);
    }
    for (i = 0; i < B->length; i++)
    {
        A->coeffs[i] = nmod_neg(B->coeffs[i], ctx->ffinfo->mod);
    }
    _nmod_mpoly_set_length(A, B->length, ctx);
}
