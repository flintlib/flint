/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod.h"
#include "n_poly.h"
#include "mpoly.h"
#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_one(fq_nmod_mpoly_t A, const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    fq_nmod_mpoly_fit_length(A, 1, ctx);
    _n_fq_one(A->coeffs + d*0, d);
    mpoly_monomial_zero(A->exps, N);
    _fq_nmod_mpoly_set_length(A, 1, ctx);
}
