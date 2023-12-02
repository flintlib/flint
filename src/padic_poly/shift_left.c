/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "padic_poly.h"

void padic_poly_shift_left(padic_poly_t rop, const padic_poly_t op, slong n,
                           const padic_ctx_t ctx)
{
    if (rop->N < op->N)
    {
        flint_throw(FLINT_ERROR, "Exception (padic_poly_shift_left).  rop->N < op->N.\n");
    }

    if (n == 0)
    {
        padic_poly_set(rop, op, ctx);
    }
    else if (op->length == 0)
    {
        padic_poly_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, op->length + n);
        _fmpz_poly_shift_left(rop->coeffs, op->coeffs, op->length, n);
        rop->val = op->val;
        _padic_poly_set_length(rop, op->length + n);
        /* TODO: Reduce */
    }
}

