/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"

/* return 1 if quotient is exact */
int fmpq_mpoly_divides(fmpq_mpoly_t Q, const fmpq_mpoly_t A,
                              const fmpq_mpoly_t B, const fmpq_mpoly_ctx_t ctx)
{
    int res;

    if (fmpq_mpoly_is_zero(B, ctx))
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in fmpq_mpoly_divides");
    }

    if (fmpq_mpoly_is_zero(A, ctx))
    {
        fmpq_mpoly_zero(Q, ctx);
        return 1;
    }

    res = fmpz_mpoly_divides(Q->zpoly, A->zpoly, B->zpoly, ctx->zctx);
    if (!res)
    {
        fmpq_mpoly_zero(Q, ctx);
        return 0;
    }

    fmpq_div(Q->content, A->content, B->content);
    return 1;
}
