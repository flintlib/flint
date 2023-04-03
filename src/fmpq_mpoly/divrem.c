/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_divrem(fmpq_mpoly_t Q, fmpq_mpoly_t R,
                          const fmpq_mpoly_t A, const fmpq_mpoly_t B,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t scale;
    fmpq_t t;

    if (fmpq_mpoly_is_zero(B, ctx))
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in fmpq_mpoly_divrem");
    }

    if (fmpq_mpoly_is_zero(A, ctx))
    {
        fmpq_mpoly_zero(Q, ctx);
        fmpq_mpoly_zero(R, ctx);
        return;
    }

    fmpz_init(scale);
    fmpz_mpoly_quasidivrem(scale, Q->zpoly, R->zpoly,
                                                A->zpoly, B->zpoly, ctx->zctx);

    fmpq_init(t);
    fmpq_div_fmpz(t, A->content, scale);
    fmpq_div(Q->content, t, B->content);
    fmpq_swap(t, R->content);
    fmpq_clear(t);
    fmpz_clear(scale);

    fmpq_mpoly_reduce(Q, ctx);
    fmpq_mpoly_reduce(R, ctx);
}
