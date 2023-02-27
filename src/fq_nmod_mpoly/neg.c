/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_neg(
    fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(B->bits, ctx->minfo);

    if (A != B)
    {
        fq_nmod_mpoly_fit_length_reset_bits(A, B->length, B->bits, ctx);
        mpoly_copy_monomials(A->exps, B->exps, B->length, N);
    }
    _nmod_vec_neg(A->coeffs, B->coeffs, d*B->length, fq_nmod_ctx_mod(ctx->fqctx));
    _fq_nmod_mpoly_set_length(A, B->length, ctx);
}
