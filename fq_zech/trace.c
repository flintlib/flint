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
fq_zech_trace(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    mp_limb_t p_i, trace;
    fq_zech_t t, op_p_i;
    double qm1inv;
    if (fq_zech_is_zero(op, ctx))
    {
        fmpz_zero(rop);
        return;
    }

    fq_zech_zero(t, ctx);
    qm1inv = n_precompute_inverse(ctx->qm1);

    for (p_i = 1; p_i <= ctx->qm1; p_i *= ctx->p)
    {
        /* op_q_i = op ^ (p ^ i) */
        op_p_i->value = n_mulmod_precomp(op->value, p_i, ctx->qm1, qm1inv);

        /* t += op_p_i */
        fq_zech_add(t, t, op_p_i, ctx);
    }

    if (fq_zech_is_zero(t, ctx))
    {
        fmpz_zero(rop);
    }
    else
    {
        trace = t->value / ctx->qm1opm1;
        trace = n_powmod(ctx->prime_root, trace, ctx->p);
        fmpz_set_ui(rop, trace);
    }
}
