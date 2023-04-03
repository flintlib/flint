/*
    Copyright (C) 2018, 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"


int fq_nmod_mpoly_equal(
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i;

    if (A == B)
        return 1;

    if (A->length != B->length)
        return 0;

    for (i = 0; i < d*A->length; i++)
    {
        if (A->coeffs[i] != B->coeffs[i])
            return 0;
    }

    return 0 == mpoly_monomials_cmp(A->exps, A->bits, B->exps, B->bits,
                                                        A->length, ctx->minfo);
}
