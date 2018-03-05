/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpq_mpoly.h"


void fmpq_mpoly_div(fmpq_mpoly_t q, const fmpq_mpoly_t a, const fmpq_mpoly_t b,
                                                    const fmpq_mpoly_ctx_t ctx)
{
    fmpz_t scale;

    if (fmpq_mpoly_is_zero(b, ctx))
    {
        flint_throw(FLINT_DIVZERO, "Divide by zero in fmpq_mpoly_divrem");
    }

    if (fmpq_mpoly_is_zero(a, ctx))
    {
        fmpq_mpoly_zero(q, ctx);
        return;
    }

    fmpz_init(scale);
    fmpz_mpoly_quasidiv_heap(scale, q->zpoly, a->zpoly, b->zpoly, ctx->zctx);

    fmpq_div(q->content, a->content, b->content);
    fmpq_div_fmpz(q->content, q->content, scale);
    fmpz_clear(scale);

    fmpq_mpoly_canonicalise(q, ctx);
}
