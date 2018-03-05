/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_divrem(fmpq_mpoly_t q, fmpq_mpoly_t r,
                  const fmpq_mpoly_t a, const fmpq_mpoly_t b,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t scale;
    fmpq_t t;

    if (fmpq_mpoly_is_zero(b, ctx))
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in fmpq_mpoly_divrem");
    }

    if (fmpq_mpoly_is_zero(a, ctx))
    {
        fmpq_mpoly_zero(q, ctx);
        fmpq_mpoly_zero(r, ctx);
        return;
    }

    fmpz_init(scale);
    fmpz_mpoly_quasidivrem_heap(scale, q->zpoly, r->zpoly,
                                                a->zpoly, b->zpoly, ctx->zctx);

    fmpq_init(t);
    fmpq_div_fmpz(t, a->content, scale);
    fmpq_div(q->content, t, b->content);
    fmpq_swap(t, r->content);
    fmpq_clear(t);
    fmpz_clear(scale);

    fmpq_mpoly_canonicalise(q, ctx);
    fmpq_mpoly_canonicalise(r, ctx);
}
