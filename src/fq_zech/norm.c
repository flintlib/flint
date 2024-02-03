/*
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fq_zech.h"

void fq_zech_norm(fmpz_t rop, const fq_zech_t op, const fq_zech_ctx_t ctx)
{
    if (fq_zech_is_zero(op, ctx))
    {
        fmpz_zero(rop);
    }
    else
    {
        fmpz_set_ui(rop, n_powmod(ctx->prime_root, op->value, ctx->p));
    }
}
