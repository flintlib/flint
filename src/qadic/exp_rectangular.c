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

void _qadic_exp_rectangular(fmpz *rop, const fmpz *op, slong v, slong len, 
                            const fmpz *a, const slong *j, slong lena, 
                            const fmpz_t p, slong N, const fmpz_t pN)
{
    const slong d = j[lena - 1];
    const slong n = _padic_exp_bound(v, N, p);

    if (n < 4)
    {
        if (n == 1)  /* y := 1 */
        {
            fmpz_one(rop);
            _fmpz_vec_zero(rop + 1, d - 1);
        }
        else if (n == 2)  /* y := 1 + x */
        {
            fmpz_t f;

            fmpz_init(f);
            fmpz_pow_ui(f, p, v);

            _fmpz_vec_scalar_mul_fmpz(rop, op, len, f);
            _fmpz_vec_zero(rop + len, d - len);
            fmpz_add_ui(rop, rop, 1);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, len, pN);

            fmpz_clear(f);
        }
        else  /* y := 1 + x + x^2/2 */
        {
            slong i;
            fmpz *x = _fmpz_vec_init(len + 1);

            fmpz_pow_ui(x + len, p, v);
            _fmpz_vec_scalar_mul_fmpz(x, op, len, x + len);

            _fmpz_poly_sqr(rop, x, len);
            if (*p != WORD(2))
            {
                for (i = 0; i < 2 * len - 1; i++)
                    if (fmpz_is_odd(rop + i))
                        fmpz_add(rop + i, rop + i, pN);
            }
            _fmpz_vec_scalar_fdiv_q_2exp(rop, rop, 2 * len - 1, 1);
            _fmpz_mod_poly_reduce(rop, 2 * len - 1, a, j, lena, pN);
            _fmpz_vec_zero(rop + (2 * len - 1), d - (2 * len - 1));
            _fmpz_mod_poly_add(rop, rop, d, x, len, pN);
            fmpz_add_ui(rop, rop, 1);
            if (fmpz_equal(rop, pN))
                fmpz_zero(rop);

            _fmpz_vec_clear(x, len + 1);
        }
    }
    else  /* n >= 4 */
    {
        const slong k = fmpz_fits_si(p) ? 
                       (n - 1 - 1) / (fmpz_get_si(p) - 1) : 0;
        const slong b = n_sqrt(n);

        slong i;
        fmpz_t c, f, pNk;
        fmpz *s, *t, *x;

        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(pNk);

        s = _fmpz_vec_init(2 * d - 1);
        t = _fmpz_vec_init(2 * d - 1);
        x = _fmpz_vec_init(d * (b + 1) + d - 1);

        fmpz_pow_ui(f, p, v);
        fmpz_pow_ui(pNk, p, N + k);

        /* Compute powers x^i of the argument */
        fmpz_one(x);
        _fmpz_vec_scalar_mul_fmpz(x + d, op, len, f);
        _fmpz_vec_zero(x + d + len, d - len);
        for (i = 2; i <= b; i++)
        {
            _fmpz_mod_poly_mul(x + i * d, x + (i - 1) * d, d, x + d, d, pNk);
            _fmpz_mod_poly_reduce(x + i * d, 2 * d - 1, a, j, lena, pNk);
        }

        _fmpz_vec_zero(rop, d);
        fmpz_one(f);

        for (i = (n + b - 1) / b - 1; i >= 0; i--)
        {
            slong lo = i * b;
            slong hi = FLINT_MIN(n - 1, lo + b - 1);

            _fmpz_vec_zero(s, d);
            fmpz_one(c);

            for ( ; hi >= lo; hi--)
            {
                _fmpz_vec_scalar_addmul_fmpz(s, x + (hi - lo) * d, d, c);
                if (hi != 0)
                    fmpz_mul_ui(c, c, hi);
            }

            _fmpz_poly_mul(t, x + b * d, d, rop, d);
            _fmpz_mod_poly_reduce(t, 2 * d - 1, a, j, lena, pNk);
            _fmpz_vec_scalar_mul_fmpz(rop, s, d, f);
            _fmpz_vec_add(rop, rop, t, d);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, pNk);

            fmpz_mul(f, f, c);
        }

        /* Note exp(x) is a unit so val(sum) == val(f). */
        i = fmpz_remove(f, f, p);
        if (i)
        {
            fmpz_pow_ui(c, p, i);
            _fmpz_vec_scalar_divexact_fmpz(rop, rop, d, c);
        }

        _padic_inv(f, f, p, N);
        _fmpz_vec_scalar_mul_fmpz(rop, rop, d, f);
        _fmpz_vec_scalar_mod_fmpz(rop, rop, d, pN);

        _fmpz_vec_clear(s, 2 * d - 1);
        _fmpz_vec_clear(t, 2 * d - 1);
        _fmpz_vec_clear(x, d * (b + 1) + d - 1);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(pNk);
    }
}

int qadic_exp_rectangular(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
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

            _qadic_exp_rectangular(t, op->coeffs, v, op->length, 
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

