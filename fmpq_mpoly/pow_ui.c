/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

void fmpq_mpoly_pow_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                           ulong k, const fmpq_mpoly_ctx_t ctx)
{
    slong kk = k;

    if (kk >= 0)
    {
        fmpq_pow_si(A->content, B->content, kk);
    }
    else if (fmpq_is_zero(B->content))
    {
        fmpq_mpoly_zero(A, ctx);
        return;
    }
    else
    {
        /* we have to work around the fact that there is no fmpq_pow_ui */
        if (!fmpq_is_pm1(B->content))
        {
            flint_throw(FLINT_ERROR, "Non-unit coefficient in fmpq_mpoly_pow_ui");
        }

        if ((k % UWORD(2)) == 0 || fmpq_is_one(B->content))
        {
            fmpq_one(A->content);
        }
        else
        {
            fmpq_one(A->content);
            fmpq_neg(A->content, A->content);
        }
    }
    fmpz_mpoly_pow_ui(A->zpoly, B->zpoly, k, ctx->zctx);
}
