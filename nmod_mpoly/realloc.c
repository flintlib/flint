/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

void nmod_mpoly_realloc(
    nmod_mpoly_t A,
    slong alloc,
    const nmod_mpoly_ctx_t ctx)
{
    slong N = mpoly_words_per_exp(A->bits, ctx->minfo);

    if (alloc == 0)             /* Clear up, reinitialise */
    {
        nmod_mpoly_clear(A, ctx);
        nmod_mpoly_init(A, ctx);
        return;
    }

    A->exps_alloc = N*alloc;
    A->exps = (ulong *) flint_realloc(A->exps, A->exps_alloc*sizeof(ulong));

    A->coeffs_alloc = alloc;
    A->coeffs = (mp_limb_t *) flint_realloc(A->coeffs, A->coeffs_alloc*sizeof(ulong));
}
