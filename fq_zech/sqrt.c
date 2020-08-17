/*
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fq_zech.h"

int
fq_zech_sqrt(fq_zech_t rop, const fq_zech_t op1, const fq_zech_ctx_t ctx)
{
    if (fq_zech_is_zero(op1, ctx) || fq_zech_is_one(op1, ctx))
    {
        rop->value = op1->value;
    } else if (ctx->p == 2)
    {
        rop->value = op1->value & 1 ? (op1->value + ctx->qm1)/2 : op1->value/2;
    } else
    {
        if (op1->value & 1) /* not a square */
            return 0;

	rop->value = op1->value/2;
    }

    return 1;
}
