/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_nmod_mpoly.h"

int fq_nmod_mpoly_cmp(
    const fq_nmod_mpoly_t A,
    const fq_nmod_mpoly_t B,
    const fq_nmod_mpoly_ctx_t ctx)
{
    slong d = fq_nmod_ctx_degree(ctx->fqctx);
    slong i;
    slong length = A->length;
    const mp_limb_t * Acoeffs = A->coeffs;
    const mp_limb_t * Bcoeffs = B->coeffs;
    int cmp;

    if (A->length != B->length)
        return A->length < B->length ? -1 : 1;

    if (length < 1)
        return 0;

    cmp = mpoly_monomials_cmp(A->exps, A->bits, B->exps, B->bits,
                                                           length, ctx->minfo);
    if (cmp != 0)
        return cmp;

    for (i = 0; i < d*length; i++)
    {
        if (Acoeffs[i] != Bcoeffs[i])
            return Acoeffs[i] < Bcoeffs[i] ? -1 : 1;
    }

    return 0;
}
