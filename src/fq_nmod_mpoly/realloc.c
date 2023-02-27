/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_realloc(
    fq_nmod_mpoly_t A,
    slong alloc,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (alloc == 0)             /* Clear up, reinitialise */
    {
        fq_nmod_mpoly_clear(A, ctx);
        fq_nmod_mpoly_init(A, ctx);
        return;
    }

    A->exps_alloc = N*alloc;
    A->exps = (ulong *) flint_realloc(A->exps, A->exps_alloc*sizeof(ulong));

    A->coeffs_alloc = d*alloc;
    A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs, A->coeffs_alloc*sizeof(ulong));
}
