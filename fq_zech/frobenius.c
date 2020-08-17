/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

void
fq_zech_frobenius(fq_zech_t rop, const fq_zech_t op, slong e,
                  const fq_zech_ctx_t ctx)
{
    double qm1_inv;
    slong d = fq_zech_ctx_degree(ctx);

    e = e % d;
    if (e < 0)
        e += d;

    if (fq_zech_is_zero(op, ctx))
    {
        fq_zech_zero(rop, ctx);
    }
    else if (e == 0)
    {
        fq_zech_set(rop, op, ctx);
    }
    else
    {
        qm1_inv = n_precompute_inverse(ctx->qm1);
        e = n_powmod(ctx->p, e, ctx->qm1);
        rop->value = n_mulmod_precomp(op->value, e, ctx->qm1, qm1_inv);
    }
}
