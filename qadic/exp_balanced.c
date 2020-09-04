/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "qadic.h"

static void 
_qadic_exp_bsplit_series(fmpz *P, fmpz_t Q, fmpz *T, 
                         const fmpz *x, slong len, slong lo, slong hi, 
                         const fmpz *a, const slong *j, slong lena)
{
    const slong d = j[lena - 1];

    if (hi - lo == 1)
    {
        _fmpz_vec_set(P, x, len);
        _fmpz_vec_zero(P + len, 2*d - 1 - len);

        fmpz_set_si(Q, lo);

        _fmpz_vec_set(T, P, 2*d - 1);
    }
    else if (hi - lo == 2)
    {
        _fmpz_poly_sqr(P, x, len);
        _fmpz_vec_zero(P + (2*len - 1), d - (2*len - 1));
        _fmpz_poly_reduce(P, 2*len - 1, a, j, lena);

        fmpz_set_si(Q, lo);
        fmpz_mul_si(Q, Q, lo + 1);

        _fmpz_vec_scalar_mul_si(T, x, len, lo + 1);
        _fmpz_vec_zero(T + len, d - len);
        _fmpz_vec_add(T, T, P, d);
    }
    else
    {
        const slong m = (lo + hi) / 2;

        fmpz *PR, *TR, *W;
        fmpz_t QR;

        PR = _fmpz_vec_init(2*d - 1);
        TR = _fmpz_vec_init(2*d - 1);
        W  = _fmpz_vec_init(2*d - 1);
        fmpz_init(QR);

        _qadic_exp_bsplit_series(P, Q, T, x, len, lo, m, a, j, lena);

        _qadic_exp_bsplit_series(PR, QR, TR, x, len, m, hi, a, j, lena);

        _fmpz_poly_mul(W, TR, d, P, d);
        _fmpz_poly_reduce(W, 2*d - 1, a, j, lena);

        _fmpz_vec_scalar_mul_fmpz(T, T, d, QR);
        _fmpz_vec_add(T, T, W, d);

        _fmpz_poly_mul(W, P, d, PR, d);
        _fmpz_poly_reduce(W, 2*d - 1, a, j, lena);
        _fmpz_vec_swap(P, W, d);

        fmpz_mul(Q, Q, QR);

        _fmpz_vec_clear(PR, 2*d - 1);
        _fmpz_vec_clear(TR, 2*d - 1);
        _fmpz_vec_clear(W,  2*d - 1);
        fmpz_clear(QR);
    }
}

static void 
_qadic_exp_bsplit(fmpz *y, const fmpz *x, slong v, slong len, 
                  const fmpz *a, const slong *j, slong lena, 
                  const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];
    const slong n = _padic_exp_bound(v, N, p);

    if (n == 1)
    {
        fmpz_one(y + 0);
        _fmpz_vec_zero(y + 1, d - 1);
    }
    else
    {
        fmpz *P, *T;
        fmpz_t Q, R;
        slong f;

        P = _fmpz_vec_init(2*d - 1);
        T = _fmpz_vec_init(2*d - 1);
        fmpz_init(Q);
        fmpz_init(R);

        _qadic_exp_bsplit_series(P, Q, T, x, len, 1, n, a, j, lena);

        fmpz_add(T + 0, T + 0, Q);  /* (T,Q) := (T,Q) + 1 */

        /* Note exp(x) is a unit so val(T) == val(Q). */
        f = fmpz_remove(Q, Q, p);
        fmpz_pow_ui(R, p, f);
        _fmpz_vec_scalar_divexact_fmpz(T, T, d, R);

        _padic_inv(Q, Q, p, N);
        _fmpz_vec_scalar_mul_fmpz(y, T, d, Q);

        _fmpz_vec_clear(P, 2*d - 1);
        _fmpz_vec_clear(T, 2*d - 1);
        fmpz_clear(Q);
        fmpz_clear(R);
    }
}

void _qadic_exp_balanced(fmpz *rop, const fmpz *x, slong v, slong len, 
                         const fmpz *a, const slong *j, slong lena, 
                         const fmpz_t p, slong N, const fmpz_t pN)
{
    const slong d = j[lena - 1];

    fmpz_t pw;
    fmpz *r, *s, *t;
    slong i, w;

    r = _fmpz_vec_init(d);
    s = _fmpz_vec_init(2*d - 1);
    t = _fmpz_vec_init(d);
    fmpz_init(pw);

    fmpz_pow_ui(pw, p, v);
    _fmpz_vec_scalar_mul_fmpz(t, x, len, pw);
    _fmpz_vec_scalar_mod_fmpz(t, t, len, pN);
    _fmpz_vec_zero(t + len, d - len);

    fmpz_set(pw, p);
    fmpz_one(rop + 0);
    _fmpz_vec_zero(rop + 1, d - 1);
    w = 1;

    while (!_fmpz_vec_is_zero(t, d))
    {
        fmpz_mul(pw, pw, pw);

        for (i = 0; i < d; i++)
        {
            fmpz_fdiv_r(r + i, t + i, pw);
            fmpz_sub(t + i, t + i, r + i);
        }

        if (!_fmpz_vec_is_zero(r, d))
        {
            _qadic_exp_bsplit(r, r, w, d, a, j, lena, p, N);
            _fmpz_poly_mul(s, rop, d, r, d);
            _fmpz_poly_reduce(s, 2*d - 1, a, j, lena);
            _fmpz_vec_scalar_mod_fmpz(rop, s, d, pN);
        }

        w *= 2;
    }

    _fmpz_vec_clear(r, d);
    _fmpz_vec_clear(s, 2*d - 1);
    _fmpz_vec_clear(t, d);
    fmpz_clear(pw);
}

int qadic_exp_balanced(qadic_t rop, const qadic_t op, const qadic_ctx_t ctx)
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

            fmpz_t pN;
            int alloc;

            alloc = _padic_ctx_pow_ui(pN, N, &ctx->pctx);

            padic_poly_fit_length(rop, d);

            _qadic_exp_balanced(rop->coeffs, op->coeffs, v, op->length, 
                                ctx->a, ctx->j, ctx->len, p, N, pN);
            rop->val = 0;

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

