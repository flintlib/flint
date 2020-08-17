/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_shift(padic_t rop, const padic_t op, slong v, const padic_ctx_t ctx)
{
    if (padic_is_zero(op) || (padic_val(op) + v >= padic_prec(rop)))
    {
        padic_zero(rop);
    }
    else
    {
        fmpz_set(padic_unit(rop), padic_unit(op));
        padic_val(rop) = padic_val(op) + v;
        _padic_reduce(rop, ctx);
    }
}

