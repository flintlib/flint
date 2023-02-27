/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic.h"

void padic_mul(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (padic_is_zero(op1) || padic_is_zero(op2) || 
        padic_val(op1) + padic_val(op2) >= padic_prec(rop))
    {
        padic_zero(rop);
    }
    else
    {
        fmpz_mul(padic_unit(rop), padic_unit(op1), padic_unit(op2));
        padic_val(rop) = padic_val(op1) + padic_val(op2);

        _padic_reduce(rop, ctx);
    }
}

