/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_set_ui(
    fq_nmod_mpoly_t A,
    ulong c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (c >= ctx->fqctx->mod.n)
        NMOD_RED(c, c, ctx->fqctx->mod);

    if (c == 0)
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    fq_nmod_mpoly_fit_length(A, 1, ctx);
    _n_fq_zero(A->coeffs + d*0, d);
    A->coeffs[0] = c;
    mpoly_monomial_zero(A->exps, N);
    _fq_nmod_mpoly_set_length(A, 1, ctx);
}


void fq_nmod_mpoly_set_fmpz(
    fq_nmod_mpoly_t A,
    const fmpz_t c,
    const fq_nmod_mpoly_ctx_t ctx)
{
    fq_nmod_mpoly_set_ui(A, fmpz_fdiv_ui(c, ctx->fqctx->mod.n), ctx);
}

