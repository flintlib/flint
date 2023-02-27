/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_poly.h"
#include "padic_poly.h"

/*
    Assumes that len > 0.
 */

static void _padic_poly_get_fmpq_poly(fmpz *rop, fmpz_t den, 
                                      const fmpz *op, slong val, slong len, 
                                      const fmpz_t p)
{
    if (val == 0)
    {
        _fmpz_vec_set(rop, op, len);
        fmpz_one(den);
    }
    else if (val == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(rop, op, len, p);
        fmpz_one(den);
    }
    else if (val > 1)
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_pow_ui(t, p, val);

        _fmpz_vec_scalar_mul_fmpz(rop, op, len, t);
        fmpz_one(den);

        fmpz_clear(t);
    }
    else
    {
        _fmpz_vec_set(rop, op, len);
        fmpz_pow_ui(den, p, -val);
    }
}

void padic_poly_get_fmpq_poly(fmpq_poly_t rop, 
                              const padic_poly_t op, const padic_ctx_t ctx)
{
    const slong len = op->length;

    if (len == 0)
    {
        fmpq_poly_zero(rop);
    }
    else
    {
        fmpq_poly_fit_length(rop, len);
        _padic_poly_get_fmpq_poly(rop->coeffs, rop->den, 
                                  op->coeffs, op->val, op->length, ctx->p);
        _fmpq_poly_set_length(rop, len);
    }
}

