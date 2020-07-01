/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_cmp(const fq_nmod_mpoly_t A, const fq_nmod_mpoly_t B,
                                                 const fq_nmod_mpoly_ctx_t ctx)
{
    int cmp;
    slong i;
    slong length = A->length;
    fq_nmod_struct * Acoeffs = A->coeffs;
    fq_nmod_struct * Bcoeffs = B->coeffs;

    if (A->length != B->length)
        return A->length < B->length ? -1 : 1;

    if (length <= 0)
        return 0;

    cmp = mpoly_monomials_cmp(A->exps, A->bits, B->exps, B->bits,
                                                           length, ctx->minfo);
    if (cmp != 0)
        return cmp;

    for (i = 0; i < length; i++)
    {
        cmp = fq_nmod_cmp(Acoeffs + i, Bcoeffs + i, ctx->fqctx);
        if (cmp != 0)
            return cmp;
    }

    return 0;
}
