/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

void fq_nmod_mpoly_realloc(fq_nmod_mpoly_t A,
                                    slong alloc, const fq_nmod_mpoly_ctx_t ctx)
{
    slong i, N;

    if (alloc == 0)             /* Clear up, reinitialise */
    {
        fq_nmod_mpoly_clear(A, ctx);
        fq_nmod_mpoly_init(A, ctx);
        return;
    }

    N = mpoly_words_per_exp(A->bits, ctx->minfo);

    for (i = alloc; i < A->alloc; i++)
        fq_nmod_clear(A->coeffs + i, ctx->fqctx);


    if (A->alloc != 0)            /* Realloc */
    {
        A->exps = (ulong *) flint_realloc(A->exps, alloc*N*sizeof(ulong));
        A->coeffs = (fq_nmod_struct *)
                        flint_realloc(A->coeffs, alloc*sizeof(fq_nmod_struct));
    }
    else                        /* Nothing allocated already so do it now */
    {
        A->exps   = (ulong *) flint_malloc(alloc*N*sizeof(ulong));
        A->coeffs = (fq_nmod_struct *)
                                    flint_malloc(alloc*sizeof(fq_nmod_struct));

    }

    for (i = A->alloc; i < alloc; i++)
        fq_nmod_init(A->coeffs + i, ctx->fqctx);

    A->alloc = alloc;
}
