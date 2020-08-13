/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_cmp(const fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    int cmp;
    slong i;
    slong length = A->zpoly->length;
    fmpz * Acoeffs = A->zpoly->coeffs;
    fmpz * Bcoeffs = B->zpoly->coeffs;

    if (A->zpoly->length != B->zpoly->length)
        return A->zpoly->length < B->zpoly->length ? -1 : 1;

    if (length <= 0)
        return 0;

    cmp = mpoly_monomials_cmp(A->zpoly->exps, A->zpoly->bits,
                              B->zpoly->exps, B->zpoly->bits,
                                                     length, ctx->zctx->minfo);
    if (cmp != 0)
        return cmp;

    cmp = fmpz_cmp(fmpq_denref(A->content), fmpq_denref(B->content));
    if (cmp != 0)
        return cmp;

    cmp = fmpz_cmp(fmpq_numref(A->content), fmpq_numref(B->content));
    if (cmp != 0)
        return cmp;

    for (i = 0; i < length; i++)
    {
        cmp = fmpz_cmp(Acoeffs + i, Bcoeffs + i);
        if (cmp != 0)
            return cmp;
    }

    return 0;
}
