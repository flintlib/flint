/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "padic.h"
#include "padic_poly.h"

void padic_poly_set_fmpq_poly(padic_poly_t rop,
                              const fmpq_poly_t op, const padic_ctx_t ctx)
{
    const slong len = op->length;

    if (len == 0)
    {
        padic_poly_zero(rop);
    }
    else
    {
        const slong N = padic_poly_prec(rop);
        fmpz_t t;

        fmpz_init(t);

        rop->val = - fmpz_remove(t, op->den, ctx->p);

        if (rop->val < N)
        {
            padic_poly_fit_length(rop, len);
            _padic_poly_set_length(rop, len);

            _padic_inv(t, t, ctx->p, N - rop->val);
            _fmpz_vec_scalar_mul_fmpz(rop->coeffs, op->coeffs, len, t);
            if (rop->val == 0)
                padic_poly_canonicalise(rop, ctx->p);

            padic_poly_reduce(rop, ctx);
        }
        else
        {
            padic_poly_zero(rop);
        }

        fmpz_clear(t);
    }
}
