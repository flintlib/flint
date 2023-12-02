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
    Computes the norm of an element $x$ of $\mathbf{Z}_q$ via the identity

        $\Norm(x) = \exp \Trace \log (x)$

    whenever $y = 1-x$ has valuation $v$ greater than $(p-1)^{-1}$.

    Assumes that $y$ is non-zero.
 */

void _qadic_norm_analytic(fmpz_t rop, const fmpz *y, slong v, slong len,
                          const fmpz *a, const slong *j, slong lena,
                          const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];
    fmpz_t pN, tru;
    slong trv;
    fmpz *lg;

    fmpz_init(pN);
    fmpz_init(tru);
    lg = _fmpz_vec_init(d);

    fmpz_pow_ui(pN, p, N);

    _qadic_log(lg, y, v, len, a, j, lena, p, N, pN);

    _qadic_trace(tru, lg, d, a, j, lena, pN);

    if (!fmpz_is_zero(tru))
    {
        trv = fmpz_remove(tru, tru, p);
        _padic_exp(rop, tru, trv, p, N);
        fmpz_mod(rop, rop, pN);
    }
    else
    {
        fmpz_one(rop);
    }

    fmpz_clear(pN);
    fmpz_clear(tru);
    _fmpz_vec_clear(lg, d);
}

void qadic_norm_analytic(padic_t rop, const qadic_t op, const qadic_ctx_t ctx)
{
    const slong N  = padic_prec(rop);
    const slong d  = qadic_ctx_degree(ctx);
    const fmpz *p = (&ctx->pctx)->p;

    /* N(p^v u) = p^{dv} N(u) */

    if (qadic_is_zero(op) || d * op->val >= N)
    {
        padic_zero(rop);
    }
    else if (op->length == 1)
    {
        fmpz_t pN;
        int alloc;

        alloc = _padic_ctx_pow_ui(pN, N - d * op->val, (&ctx->pctx));

        fmpz_powm_ui(padic_unit(rop), op->coeffs + 0, d, pN);
        padic_val(rop) = d * op->val;

        if (alloc)
            fmpz_clear(pN);
    }
    else  /* len >= 2 */
    {
        fmpz *y;
        slong w;

        y = _fmpz_vec_init(op->length);

        _fmpz_vec_neg(y, op->coeffs, op->length);
        fmpz_add_ui(y + 0, y + 0, 1);
        w = _fmpz_vec_ord_p(y, op->length, p);

        if ((w < 2 && *p == WORD(2)) || w < 1)
        {
            flint_throw(FLINT_ERROR, "ERROR (qadic_norm_analytic).  w = %wd.\n", w);
        }

        _qadic_norm_analytic(padic_unit(rop), y, w, op->length,
                             ctx->a, ctx->j, ctx->len, p, N - d * op->val);
        padic_val(rop) = d * op->val;

        _fmpz_vec_clear(y, op->length);
    }
}

