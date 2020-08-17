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
fq_zech_mul_si(fq_zech_t rop, const fq_zech_t op, slong x,
               const fq_zech_ctx_t ctx)
{
    mp_limb_t y;
    if (x == 0 || fq_zech_is_zero(op, ctx))
    {
        fq_zech_zero(rop, ctx);
        return;
    }
    if (x < 0)
    {
        y = -x;
        fq_zech_mul_ui(rop, op, y, ctx);
        fq_zech_neg(rop, rop, ctx);
    }
    else
    {
        y = x;
        fq_zech_mul_ui(rop, op, y, ctx);
    }
}
