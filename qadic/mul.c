/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qadic.h"

/*
    Forms the product of (op1,len1) and (op2,len2) modulo (a,j,lena) and pN.
    Assumes that len1 >= len2 > 0.  Requires rop to be of size at least 
    len1 + len2 - 1.
 */

static 
void _qadic_mul(fmpz *rop, 
                const fmpz *op1, slong len1, const fmpz *op2, slong len2, 
                const fmpz *a, const slong *j, slong lena, const fmpz_t pN)
{
    _fmpz_poly_mul(rop, op1, len1, op2, len2);
    _fmpz_mod_poly_reduce(rop, len1 + len2 - 1, a, j, lena, pN);
}

void qadic_mul(qadic_t x, const qadic_t y, const qadic_t z, 
                          const qadic_ctx_t ctx)
{
    const slong leny = y->length;
    const slong lenz = z->length;
    const slong lenx = leny + lenz - 1;
    const slong N    = qadic_prec(x);
    const slong d    = qadic_ctx_degree(ctx);

    if (leny == 0 || lenz == 0 || y->val + z->val >= N)
    {
        qadic_zero(x);
    }
    else
    {
        fmpz *t;
        fmpz_t pN;
        int alloc;

        x->val = y->val + z->val;

        alloc = _padic_ctx_pow_ui(pN, N - x->val, &ctx->pctx);

        if (x == y || x == z)
        {
            t = _fmpz_vec_init(lenx);
        }
        else
        {
            padic_poly_fit_length(x, lenx);
            t = x->coeffs;
        }

        if (leny >= lenz)
            _qadic_mul(t, y->coeffs, leny, 
                          z->coeffs, lenz, ctx->a, ctx->j, ctx->len, pN);
        else
            _qadic_mul(t, z->coeffs, lenz, 
                          y->coeffs, leny, ctx->a, ctx->j, ctx->len, pN);

        if (x == y || x == z)
        {
            _fmpz_vec_clear(x->coeffs, x->alloc);
            x->coeffs = t;
            x->alloc  = lenx;
        }

        _padic_poly_set_length(x, FLINT_MIN(lenx, d));
        _padic_poly_normalise(x);

        if (alloc)
            fmpz_clear(pN);
    }
}

