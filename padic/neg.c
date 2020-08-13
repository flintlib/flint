/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_neg(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op) || padic_val(op) >= padic_prec(rop))
    {
        padic_zero(rop);
    }
    else
    {
        fmpz_t pow;
        int alloc;

        padic_val(rop) = padic_val(op);

        alloc = _padic_ctx_pow_ui(pow, padic_prec(rop) - padic_val(rop), ctx);
        fmpz_sub(padic_unit(rop), pow, padic_unit(op));
        if (alloc)
            fmpz_clear(pow);

        _padic_reduce(rop, ctx);
    }
}

