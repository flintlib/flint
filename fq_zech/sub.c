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
fq_zech_sub(fq_zech_t rop, const fq_zech_t op1, const fq_zech_t op2,
            const fq_zech_ctx_t ctx)
{
    mp_limb_t index, c;
    if (op2->value == ctx->qm1)
    {
        rop->value = op1->value;
    }
    else if (op1->value == ctx->qm1)
    {
        fq_zech_neg(rop, op2, ctx);
    }
    else
    {
        index = n_submod(op2->value, op1->value, ctx->qm1);
        index = n_submod(index, ctx->qm1o2, ctx->qm1);

        c = ctx->zech_log_table[index];
        if (c != ctx->qm1)
        {
            c = n_addmod(c, op1->value, ctx->qm1);
        }
        rop->value = c;
    }
}
