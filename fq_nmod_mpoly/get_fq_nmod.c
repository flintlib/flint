/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_get_fq_nmod(fq_nmod_t c, const fq_nmod_mpoly_t A, 
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    slong N;

    if (A->length > WORD(1))
        flint_throw(FLINT_ERROR, "Nonconstant polynomial in fq_nmod_mpoly_get_fq_nmod");

    if (A->length == WORD(0))
    {
        fq_nmod_zero(c, ctx->fqctx);
        return;
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);
    if (!mpoly_monomial_is_zero(A->exps + N*0, N))
        flint_throw(FLINT_ERROR, "Nonconstant monomial in fq_nmod_mpoly_get_fq_nmod");

    fq_nmod_set(c, A->coeffs + 0, ctx->fqctx);
}
