/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_set_ui(fq_nmod_mpoly_t A, ulong c,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (c >= ctx->fqctx->modulus->mod.n)
    {
        NMOD_RED(c, c, ctx->fqctx->modulus->mod);
    }

    if (c == UWORD(0))
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    fq_nmod_mpoly_fit_length(A, 1, ctx);
    fq_nmod_set_ui(A->coeffs + 0, c, ctx->fqctx);
    FLINT_ASSERT(!fq_nmod_is_zero(A->coeffs, ctx->fqctx));
    mpoly_monomial_zero(A->exps, N);
    _fq_nmod_mpoly_set_length(A, 1, ctx);
}
