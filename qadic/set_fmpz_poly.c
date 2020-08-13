/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qadic.h"

void qadic_set_fmpz_poly(qadic_t rop, const fmpz_poly_t op, 
                         const qadic_ctx_t ctx)
{
    const slong len = op->length;

    if (len == 0)
    {
        qadic_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, len);
        _padic_poly_set_length(rop, len);
        _fmpz_vec_set(rop->coeffs, op->coeffs, len);
        rop->val = 0;

        qadic_reduce(rop, ctx);
    }
}

