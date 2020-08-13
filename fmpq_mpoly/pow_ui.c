/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

int fmpq_mpoly_pow_ui(fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                           ulong k, const fmpq_mpoly_ctx_t ctx)
{
    slong success;

    /* we have to work around the fact that there is no fmpq_pow_ui */
    if (k <= WORD_MAX)
    {
        fmpq_pow_si(A->content, B->content, k);
        success = 1;
    }
    else if (fmpq_is_zero(B->content))
    {
        fmpq_mpoly_zero(A, ctx);
        return 1;
    }
    else
    {
        if (!fmpq_is_pm1(B->content))
        {
            success = 0;
        }
        else
        {
            fmpz_set_si(fmpq_numref(A->content),
                      (k % UWORD(2)) == 0 || fmpq_is_one(B->content) ? 1 : -1);
            fmpz_one(fmpq_denref(A->content));
            success = 1;
        }
    }

    success = success && fmpz_mpoly_pow_ui(A->zpoly, B->zpoly, k, ctx->zctx);

    return success;
}
