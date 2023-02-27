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
fq_zech_pth_root(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
{
    slong i, d;
    mp_limb_t e;
    double qm1inv;

    if (fq_zech_is_zero(op1, ctx) || fq_zech_is_one(op1, ctx))
    {
        rop->value = op1->value;
        return;
    }
    
    d = fq_zech_ctx_degree(ctx) - 1;
    qm1inv = n_precompute_inverse(ctx->qm1);

    e = op1->value;
    for (i = 0; i < d; i++)
    {
        e = n_mulmod_precomp(ctx->p, e, ctx->qm1, qm1inv);
    }
    
    rop->value = e;    
}
