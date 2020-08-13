/*
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpq_poly.h"
#include "padic_poly.h"

void padic_poly_set_fmpq_poly(padic_poly_t f, 
                              const fmpq_poly_t g, const padic_ctx_t ctx)
{
    const slong len = g->length;

    if (len == 0)
    {
        padic_poly_zero(f);
    }
    else
    {
        const slong N = padic_poly_prec(f);
        fmpz_t t;

        fmpz_init(t);

        f->val = - fmpz_remove(t, g->den, ctx->p);

        if (f->val < N)
        {
            padic_poly_fit_length(f, len);
            _padic_poly_set_length(f, len);

            _padic_inv(t, t, ctx->p, N - f->val);
            _fmpz_vec_scalar_mul_fmpz(f->coeffs, g->coeffs, len, t);
            if (f->val == 0)
                padic_poly_canonicalise(f, ctx->p);

            padic_poly_reduce(f, ctx);
        }
        else
        {
            padic_poly_zero(f);
        }

        fmpz_clear(t);
    }
}

