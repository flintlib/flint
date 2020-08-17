/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
_fq_trace(fmpz_t rop, const fmpz * op, slong len, const fq_ctx_t ctx)
{
    const slong d = fq_ctx_degree(ctx);

    slong i, l;
    fmpz *t;

    t = _fmpz_vec_init(d);

    fmpz_set_ui(t + 0, d);
    for (i = 1; i < d; i++)
    {
        for (l = ctx->len - 2; l >= 0 && ctx->j[l] >= d - (i - 1); l--)
        {
            fmpz_addmul(t + i, t + (ctx->j[l] + i - d), ctx->a + l);
        }

        if (l >= 0 && ctx->j[l] == d - i)
        {
            fmpz_addmul_ui(t + i, ctx->a + l, i);
        }

        fmpz_neg(t + i, t + i);
        fmpz_mod(t + i, t + i, fq_ctx_prime(ctx));
    }

    fmpz_zero(rop);
    for (i = 0; i < len; i++)
    {
        fmpz_addmul(rop, op + i, t + i);
    }
    fmpz_mod(rop, rop, fq_ctx_prime(ctx));

    _fmpz_vec_clear(t, d);
}

void
fq_trace(fmpz_t rop, const fq_t op, const fq_ctx_t ctx)
{
    if (fq_is_zero(op, ctx))
    {
        fmpz_zero(rop);
    }
    else
    {
        _fq_trace(rop, op->coeffs, op->length, ctx);
    }
}
