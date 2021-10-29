/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly.h"

int fmpz_mod_mpoly_cmp(
    const fmpz_mod_mpoly_t A,
    const fmpz_mod_mpoly_t B,
    const fmpz_mod_mpoly_ctx_t ctx)
{
    int cmp;
    slong i;
    slong length = A->length;
    fmpz * Acoeffs = A->coeffs;
    fmpz * Bcoeffs = B->coeffs;

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
        cmp = fmpz_cmp(Acoeffs + i, Bcoeffs + i);
        if (cmp != 0)
            return cmp < 0 ? -1 : 1;
    }

    return 0;
}
