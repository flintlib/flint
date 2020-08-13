/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_set_fq_nmod(fq_nmod_mpoly_t A,
                              const fq_nmod_t c, const fq_nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (fq_nmod_is_zero(c, ctx->fqctx))
    {
        fq_nmod_mpoly_zero(A, ctx);
        return;
    }

    fq_nmod_mpoly_fit_length(A, 1, ctx);
    fq_nmod_set(A->coeffs + 0, c, ctx->fqctx);
    mpoly_monomial_zero(A->exps + N*0, N);
    _fq_nmod_mpoly_set_length(A, 1, ctx);
}
