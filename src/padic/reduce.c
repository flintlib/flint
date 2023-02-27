/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void _padic_reduce(padic_t rop, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(padic_unit(rop)))
    {
        if (padic_val(rop) < padic_prec(rop))
        {
            int c;
            fmpz_t pow;

            c = _padic_ctx_pow_ui(pow, padic_prec(rop) - padic_val(rop), ctx);
            fmpz_mod(padic_unit(rop), padic_unit(rop), pow);
            if (c)
                fmpz_clear(pow);
        }
        else
        {
            padic_zero(rop);
        }
    }
}

void padic_reduce(padic_t rop, const padic_ctx_t ctx)
{
    _padic_canonicalise(rop, ctx);
    _padic_reduce(rop, ctx);
}

