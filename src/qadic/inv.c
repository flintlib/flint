/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "qadic.h"

static int __fmpz_mod_poly_invmod(fmpz *A,
                          const fmpz *B, slong lenB,
                          const fmpz *P, slong lenP, const fmpz_t p)
{
    fmpz * t, * u;
    fmpz_mod_ctx_t mod;
    int res;

    t = _fmpz_vec_init(lenB);
    u = _fmpz_vec_init(lenP);

    fmpz_mod_ctx_init(mod, p);

    _fmpz_vec_scalar_mod_fmpz(t, B, lenB, p);
    _fmpz_vec_scalar_mod_fmpz(u, P, lenP, p);

    res = _fmpz_mod_poly_invmod(A, t, lenB, u, lenP, mod);

    fmpz_mod_ctx_clear(mod);

    _fmpz_vec_clear(t, lenB);
    _fmpz_vec_clear(u, lenP);

    return res;
}

void _qadic_inv(fmpz *rop, const fmpz *op, slong len,
                const fmpz *a, const slong *j, slong lena,
                const fmpz_t p, slong N)
{
    const slong d = j[lena - 1];

    if (len == 1)
    {
        _padic_inv(rop, op, p, N);
        _fmpz_vec_zero(rop + 1, d - 1);
    }
    else if (N == 1)
    {
        fmpz *P = _fmpz_vec_init(d + 1);
        slong k;

        for (k = 0; k < lena; k++)
            fmpz_set(P + j[k], a + k);

        __fmpz_mod_poly_invmod(rop, op, len, P, d + 1, p);

        _fmpz_vec_clear(P, d + 1);
    }
    else  /* d, N >= 2 */
    {
        slong *e, i, n;
        fmpz *pow, *u;
        fmpz *s, *t;

        n = FLINT_CLOG2(N) + 1;

        /* Compute sequence of exponents */
        e = flint_malloc(n * sizeof(slong));
        for (e[i = 0] = N; e[i] > 1; i++)
            e[i + 1] = (e[i] + 1) / 2;

        pow = _fmpz_vec_init(n);
        u   = _fmpz_vec_init(len * n);
        s   = _fmpz_vec_init(2 * d - 1);
        t   = _fmpz_vec_init(2 * d - 1);

        /* Compute powers of p */
        {
            fmpz_one(t);
            fmpz_set(pow + i, p);
        }
        for (i--; i >= 1; i--)
        {
            if (e[i] & WORD(1))
            {
                fmpz_mul(pow + i, t, pow + (i + 1));
                fmpz_mul(t, t, t);
            }
            else
            {
                fmpz_mul(t, t, pow + (i + 1));
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
            }
        }
        {
            if (e[i] & WORD(1))
                fmpz_mul(pow + i, t, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }

        /* Compute reduced units */
        {
            _fmpz_vec_scalar_mod_fmpz(u + 0 * len, op, len, pow + 0);
        }
        for (i = 1; i < n; i++)
        {
            _fmpz_vec_scalar_mod_fmpz(u + i * len, u + (i - 1) * len, len, pow + i);
        }

        /* Run Newton iteration */
        i = n - 1;
        {
            fmpz *P = _fmpz_vec_init(d + 1);
            slong k;

            for (k = 0; k < lena; k++)
                fmpz_set(P + j[k], a + k);

            __fmpz_mod_poly_invmod(rop, u + i * len, len, P, d + 1, pow + i);

            _fmpz_vec_clear(P, d + 1);
        }
        for (i--; i >= 0; i--)  /* z' := 2 z - a z^2 */
        {
            _fmpz_poly_sqr(s, rop, d);
            _fmpz_poly_reduce(s, 2 * d - 1, a, j, lena);

            _fmpz_poly_mul(t, s, d, u + i * len, len);
            _fmpz_poly_reduce(t, d + len - 1, a, j, lena);

            _fmpz_vec_scalar_mul_2exp(rop, rop, d, 1);
            _fmpz_poly_sub(rop, rop, d, t, d);
            _fmpz_vec_scalar_mod_fmpz(rop, rop, d, pow + i);
        }

        _fmpz_vec_clear(pow, n);
        _fmpz_vec_clear(u, len * n);
        _fmpz_vec_clear(s, 2 * d - 1);
        _fmpz_vec_clear(t, 2 * d - 1);
        flint_free(e);
    }
}

void qadic_inv(qadic_t x, const qadic_t y, const qadic_ctx_t ctx)
{
    const slong N = qadic_prec(x);

    if (qadic_is_zero(y))
    {
        flint_throw(FLINT_ERROR, "Exception (qadic_inv).  Zero is not invertible.\n");
    }

    /*
        If y = u p^v has negative valuation with N <= -v then the
        exact inverse of y is zero when reduced modulo $p^N$
     */
    if (N + y->val <= 0)
    {
        qadic_zero(x);
    }
    else
    {
        const slong d = qadic_ctx_degree(ctx);
        fmpz *t;

        if (x == y)
        {
            t = _fmpz_vec_init(d);
        }
        else
        {
            padic_poly_fit_length(x, d);
            t = x->coeffs;
        }

        _qadic_inv(t, y->coeffs, y->length,
                   ctx->a, ctx->j, ctx->len, (&ctx->pctx)->p, N + y->val);
        x->val = - y->val;

        if (x == y)
        {
            _fmpz_vec_clear(x->coeffs, x->alloc);
            x->coeffs = t;
            x->alloc  = d;
            x->length = d;
        }
        else
        {
            _padic_poly_set_length(x, d);
        }
        _padic_poly_normalise(x);
    }
}

