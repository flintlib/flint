/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "padic_poly.h"

void _padic_poly_derivative(fmpz *rop, slong *rval, slong N, 
                            const fmpz *op, slong val, slong len, 
                            const padic_ctx_t ctx)
{
    fmpz_t pow;
    int alloc;

    _fmpz_poly_derivative(rop, op, len);
    *rval = val;

    alloc = _padic_ctx_pow_ui(pow, N - *rval, ctx);

    _fmpz_vec_scalar_mod_fmpz(rop, rop, len - 1, pow);
    _padic_poly_canonicalise(rop, rval, len - 1, ctx->p);

    if (alloc)
        fmpz_clear(pow);
}

void padic_poly_derivative(padic_poly_t rop, 
                           const padic_poly_t op, const padic_ctx_t ctx)
{
    const slong len = op->length;

    if (len < 2 || op->val >= rop->N)
    {
        padic_poly_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, len - 1);
        _padic_poly_derivative(rop->coeffs, &(rop->val), rop->N, 
                               op->coeffs, op->val, len, ctx);
        _padic_poly_set_length(rop, len - 1);
        _padic_poly_normalise(rop);
    }
}

