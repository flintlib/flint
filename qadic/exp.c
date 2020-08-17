/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_poly.h"
#include "qadic.h"

void _qadic_exp(fmpz *rop, const fmpz *op, slong v, slong len, 
                           const fmpz *a, const slong *j, slong lena, 
                           const fmpz_t p, slong N, const fmpz_t pN)
{
    if (N < (WORD(1) << 13) / (slong) fmpz_bits(p))
    {
        _qadic_exp_rectangular(rop, op, v, len, a, j, lena, p, N, pN);
    }
    else
    {
        const slong d = j[lena - 1];

        _qadic_exp_balanced(rop, op, v, len, a, j, lena, p, N, pN);
        _fmpz_vec_zero(rop + d, d - 1);
    }
}

int qadic_exp(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const slong N  = qadic_prec(rop);
    const slong v  = op->val;
    const fmpz *p = (&ctx->pctx)->p;

    if (padic_poly_is_zero(op))
    {
        padic_poly_one(rop);
        return 1;
    }

    if ((*p == WORD(2) && v <= 1) || (v <= 0))
    {
        return 0;
    }
    else
    {
        if (v < N)
        {
            const slong d = qadic_ctx_degree(ctx);

            fmpz *t;
            fmpz_t pN;
            int alloc;

            alloc = _padic_ctx_pow_ui(pN, N, &ctx->pctx);

            if (rop == op)
            {
                t = _fmpz_vec_init(2 * d - 1);
            }
            else
            {
                padic_poly_fit_length(rop, 2 * d - 1);
                t = rop->coeffs;
            }

            _qadic_exp(t, op->coeffs, v, op->length, 
                       ctx->a, ctx->j, ctx->len, p, N, pN);
            rop->val = 0;

            if (rop == op)
            {
                _fmpz_vec_clear(rop->coeffs, rop->alloc);
                rop->coeffs = t;
                rop->alloc  = 2 * d - 1;
                rop->length = d;
            }
            _padic_poly_set_length(rop, d);
            _padic_poly_normalise(rop);

            if (alloc)
                fmpz_clear(pN);
        }
        else
        {
            padic_poly_one(rop);
        }
        return 1;
    }
}

