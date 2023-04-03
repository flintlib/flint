/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

int nmod_mpoly_cmp(const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int cmp;
    slong i;
    slong length = A->length;
    mp_limb_t * Acoeffs = A->coeffs;
    mp_limb_t * Bcoeffs = B->coeffs;

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
        if (Acoeffs[i] != Bcoeffs[i])
            return Acoeffs[i] < Bcoeffs[i] ? -1 : 1;
    }

    return 0;
}
