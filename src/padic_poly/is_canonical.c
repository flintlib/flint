/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

int _padic_poly_is_canonical(const fmpz *op, slong val, slong len, 
                             const padic_ctx_t ctx)
{
    if (len == 0)
    {
        return (val == 0);
    }
    else
    {
        slong w = _fmpz_vec_ord_p(op, len, ctx->p);

        return (w == 0);
    }
}

int padic_poly_is_canonical(const padic_poly_t op, const padic_ctx_t ctx)
{
    return _padic_poly_is_canonical(op->coeffs, op->val, op->length, ctx);
}

