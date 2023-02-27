/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>

#include "fq_zech.h"

void
fq_zech_set_fq_nmod(fq_zech_t rop, const fq_nmod_t op, const fq_zech_ctx_t ctx)
{
    mp_limb_t i;
    fq_zech_t t;
    fq_zech_zero(rop, ctx);
    for (i = 0; i < op->length; i++)
    {
        if (op->coeffs[i] == 0)
        {
            continue;
        }
        t->value = i;
        fq_zech_mul_ui(t, t, op->coeffs[i], ctx);
        fq_zech_add(rop, rop, t, ctx);
    }
}
