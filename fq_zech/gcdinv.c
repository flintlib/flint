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
fq_zech_gcdinv(fq_zech_t rop, fq_zech_t inv, const fq_zech_t op,
               const fq_zech_ctx_t ctx)
{
    if (fq_zech_is_zero(op, ctx))
    {
        fq_zech_zero(rop, ctx);
        return;
    }
    fq_zech_one(rop, ctx);
    fq_zech_inv(inv, op, ctx);
}
