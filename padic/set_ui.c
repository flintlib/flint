/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_set_ui(padic_t rop, ulong op, const padic_ctx_t ctx)
{
    if (op == 0)
    {
        padic_zero(rop);
    }
    else if (fmpz_cmp_ui(ctx->p, op) > 0)
    {
        fmpz_set_ui(padic_unit(rop), op);
        padic_val(rop) = 0;
    }
    else
    {
        ulong p = fmpz_get_ui(ctx->p), q, r;

        /* Remove factors of p */
        padic_val(rop) = 0;
        r = n_divrem2_precomp(&q, op, p, ctx->pinv);
        while (r == 0)
        {
            op = q;
            padic_val(rop)++;
            r = n_divrem2_precomp(&q, op, p, ctx->pinv);
        }

        fmpz_set_ui(padic_unit(rop), op);

        _padic_reduce(rop, ctx);
    }
}

