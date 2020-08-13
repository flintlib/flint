/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq.h"

void
fq_pth_root(fq_t rop, const fq_t op1, const fq_ctx_t ctx)
{
    slong i, d;
    if (fq_is_zero(op1, ctx) || fq_is_one(op1, ctx))
    {
        fq_set(rop, op1, ctx);
        return;
    }

    d = fq_ctx_degree(ctx) - 1;
    fq_set(rop, op1, ctx);
    for (i = 0; i < d; i++)
    {
        fq_pow(rop, rop, fq_ctx_prime(ctx), ctx);
    }
}
