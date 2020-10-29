/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_get_term_monomial(
    fq_nmod_mpoly_t M,
    const fq_nmod_mpoly_t A,
    slong i,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    flint_bitcnt_t bits = A->bits;
    slong N = mpoly_words_per_exp(bits, ctx->minfo);

    if (i >= (ulong) A->length)
    {
        flint_throw(FLINT_ERROR, "fq_nmod_mpoly_get_term_monomial: index out of range");
    }

    fq_nmod_mpoly_fit_length_reset_bits(M, 1, bits, ctx);

    mpoly_monomial_set(M->exps + N*0, A->exps + N*i, N);
    _n_fq_one(M->coeffs + d*0, d);
    _fq_nmod_mpoly_set_length(M, 1, ctx);
}
