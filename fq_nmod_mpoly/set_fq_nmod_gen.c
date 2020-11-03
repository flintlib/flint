/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_set_fq_nmod_gen(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(A, 1, ctx);
    n_fq_gen(A->coeffs + 0, ctx->fqctx);
    mpoly_monomial_zero(A->exps, N);
    _fq_nmod_mpoly_set_length(A, !_n_fq_is_zero(A->coeffs + 0, d), ctx);
}
