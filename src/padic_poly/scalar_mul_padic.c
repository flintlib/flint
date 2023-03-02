/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

void _padic_poly_scalar_mul_padic(fmpz *rop, slong *rval, slong N, 
                                  const fmpz *op, slong val, slong len, 
                                  const padic_t c, const padic_ctx_t ctx)
{
    if (padic_is_zero(c) || val + padic_val(c) >= N)
    {
        _fmpz_vec_zero(rop, len);
        *rval = 0;
    }
    else
    {
        fmpz_t pow;
        int alloc;

        *rval = val + padic_val(c);

        alloc = _padic_ctx_pow_ui(pow, N - *rval, ctx);

        _fmpz_vec_scalar_mul_fmpz(rop, op, len, padic_unit(c));

        _fmpz_vec_scalar_mod_fmpz(rop, rop, len, pow);

        if (alloc)
            fmpz_clear(pow);    
    }
}

void padic_poly_scalar_mul_padic(padic_poly_t rop, const padic_poly_t op, 
                                 const padic_t c, const padic_ctx_t ctx)
{
    if (padic_poly_is_zero(op) || padic_is_zero(c) ||
        op->val + padic_val(c) >= rop->N)
    {
        padic_poly_zero(rop);
    }
    else
    {
        padic_poly_fit_length(rop, op->length);
        _padic_poly_set_length(rop, op->length);

        _padic_poly_scalar_mul_padic(rop->coeffs, &(rop->val), rop->N, 
                                     op->coeffs, op->val, op->length, c, ctx);
    }
}

