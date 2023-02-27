/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "padic_poly.h"

int padic_poly_get_fmpz_poly(fmpz_poly_t rop, const padic_poly_t op, 
                             const padic_ctx_t ctx)
{
    const slong len = op->length;

    if (op->val < 0)
    {
        return 0;
    }

    if (padic_poly_is_zero(op))
    {
        fmpz_poly_zero(rop);
        return 1;
    }

    fmpz_poly_fit_length(rop, len);
    _fmpz_poly_set_length(rop, len);

    if (op->val == 0)
    {
        _fmpz_vec_set(rop->coeffs, op->coeffs, len);
    }
    else  /* op->val > 0 */
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, op->val);
        _fmpz_vec_scalar_mul_fmpz(rop->coeffs, op->coeffs, len, pow);
        fmpz_clear(pow);
    }

    return 1;
}

