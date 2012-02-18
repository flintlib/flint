/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "padic.h"

/*
    Returns an integer $i$ such that for all $j \geq i$ with $j$ a 
    word-sized integer, we have $\ord_p(x^j / j!) \geq N$, where 
    $\ord_p(x) = v$.

    When $p$ is a word-sized prime, 
    returns $\ceil{\frac{(p-1)N - 1}{(p-1)v - 1}}$.
    Otherwise, returns $\ceil{N/v}$.

    Assumes that $v < N$.
 */
static long _padic_exp_bound(long v, long N, const fmpz_t p)
{
    if (fmpz_fits_si(p))
    {
        fmpz_t n, d, f;
        long i;

        fmpz_init(n);
        fmpz_init(d);
        fmpz_init(f);

        fmpz_sub_ui(f, p, 1);
        fmpz_mul_ui(n, f, N);
        fmpz_sub_ui(n, n, 1);
        fmpz_mul_ui(d, f, v);
        fmpz_sub_ui(d, d, 1); 

        fmpz_cdiv_q(f, n, d);
        i = fmpz_get_si(f);

        fmpz_clear(n);
        fmpz_clear(d);
        fmpz_clear(f);

        return i;
    }
    else
    {
        return (N + (v - 1)) / v;
    }
}

static void
padic_exp_bsplit_series(fmpz_t P, fmpz_t Q, fmpz_t T,
                        const fmpz_t x, long a, long b)
{
    if (b - a == 1)
    {
        fmpz_set(P, x);
        fmpz_set_ui(Q, a);
        fmpz_set(T, x);
    }
    else if (b - a == 2)
    {
        fmpz_mul(P, x, x);
        fmpz_set_ui(Q, a);
        fmpz_mul_ui(Q, Q, a + 1);
        fmpz_mul_ui(T, x, a + 1);
        fmpz_add(T, T, P);
    }
    else
    {
        const long m = (a + b) / 2;

        fmpz_t PR, QR, TR;

        fmpz_init(PR);
        fmpz_init(QR);
        fmpz_init(TR);

        padic_exp_bsplit_series(P, Q, T, x, a, m);
        padic_exp_bsplit_series(PR, QR, TR, x, m, b);

        fmpz_mul(T, T, QR);
        fmpz_addmul(T, P, TR);
        fmpz_mul(P, P, PR);
        fmpz_mul(Q, Q, QR);

        fmpz_clear(PR);
        fmpz_clear(QR);
        fmpz_clear(TR);
    }
}

/*
    Assumes that exponential series converges at $x$.  In particular, 
    this implies that $\ord_p(x) \geq 1$.
 */

void
padic_exp_bsplit(padic_t y, const padic_t x, const padic_ctx_t ctx)
{
    const long n = _padic_exp_bound(padic_val(x), ctx->N, ctx->p);

    if (n == 1)
    {
        padic_one(y, ctx);
    }
    else
    {
        fmpz_t t, P, Q, T;

        fmpz_init(t);
        fmpz_init(P);
        fmpz_init(Q);
        fmpz_init(T);

        padic_get_fmpz(t, x, ctx);

        padic_exp_bsplit_series(P, Q, T, t, 1, n);

        fmpz_add(T, T, Q);  /* (T,Q) := (T,Q) + 1 */

        /* As exp(x) is a unit, ord_p(T) == ord_p(Q) */
        if (_fmpz_remove(T, ctx->p, ctx->pinv))
            _fmpz_remove(Q, ctx->p, ctx->pinv);

        _padic_inv(Q, Q, ctx->p, ctx->N);  /* XXX: Aliasing */
        fmpz_mul(padic_unit(y), T, Q);
        padic_val(y) = 0;
        _padic_reduce(y, ctx);

        fmpz_clear(t);
        fmpz_clear(P);
        fmpz_clear(Q);
        fmpz_clear(T);
    }
}

/*
    XXX:  Assumes that $\ord_p(x) = 1$.
 */
void
padic_exp_balanced(padic_t y, const padic_t x, const padic_ctx_t ctx)
{
    fmpz_t t, ppow;
    padic_t r;
    long v, w;

    fmpz_init_set(ppow, ctx->p);
    fmpz_init_set(t, padic_unit(x));
    padic_init(r, ctx);

    padic_one(y, ctx);

    for (v = 2, w = 1; v < (2 * ctx->N); v *= 2)
    {
        fmpz_mul(ppow, ppow, ppow);
        fmpz_fdiv_qr(t, padic_unit(r), t, ppow);
        padic_val(r) = w;
        padic_exp_bsplit(r, r, ctx);
        padic_mul(y, y, r, ctx);
        w += v;
    }

    fmpz_clear(ppow);
    fmpz_clear(t);
    padic_clear(r, ctx);
}

/*
    TODO:  This implementation is currently broken and does not pass
    the test code for the exponential function.
 */

void
padic_exp_rectangular(padic_t res, const padic_t x, const padic_ctx_t ctx)
{
    long i, j, u, m, workprec, nterms, nxpow, nsums;

    fmpz * xpow;
    fmpz_t sum, c, f, s, t;
    fmpz_t mod;
    int alloc;

    /* XXX: Check precisions */
    i = _padic_exp_bound(padic_val(x), ctx->N, ctx->p);
    j = fmpz_fits_si(ctx->p) ? (i - 1) / (fmpz_get_si(ctx->p) - 1) : 0;
    nterms = i + 1;
    workprec = ctx->N + j;

    alloc = _padic_ctx_pow_ui(mod, workprec, ctx);

    nxpow = n_sqrt(nterms);
    nsums = (nterms + (nxpow - 1)) / nxpow;
    u     = nterms - 1;

    /*
        printf("N=%ld nterms=%d workprec=%d nxpow=%d nsums=%d\n", ctx->N,
                nterms, workprec, nxpow, nsums);
     */

    fmpz_init(sum);
    fmpz_init(c);
    fmpz_init(f);
    fmpz_init(s);
    fmpz_init(t);

    /* Compute powers: xpow[i] == x^i */
    xpow = _fmpz_vec_init(nxpow);
    padic_get_fmpz(xpow + 1, x, ctx);
    for (i = 2; i < nxpow; i++)
    {
        fmpz_mul(xpow + i, xpow + i - 1, xpow + 1);
        fmpz_mod(xpow + i, xpow + i, mod);
    }

    /* Compute x^nxpow */
    fmpz_mul(xpow, xpow + nxpow - 1, xpow + 1);
    fmpz_mod(xpow, xpow, mod);

    for (i = nsums - 1; i >= 0; i--)
    {
        /*
            Special case: only 1 term in last row;
            TODO: Handle more cleanly
         */
        if (i * nxpow == nterms - 1)
        {
            fmpz_set(sum, xpow);
            fmpz_set_ui(f, nterms - 1);
            u--;
        }
        else
        {
            const long bot = i * nxpow + 1;

            fmpz_set(s, xpow + nxpow - 1);
            fmpz_set_ui(c, u);

            while (u > bot)
            {
                fmpz_addmul(s, xpow + u - bot, c);
                u--;
                fmpz_mul_ui(c, c, u);
            }

            u -= 2;
            fmpz_add(s, s, c);

            if (i == nsums - 1)
            {
                fmpz_set(f, c);
                fmpz_set(sum, s);
            }
            else
            {
                fmpz_mul_ui(f, f, (i + 1) * nxpow);
                fmpz_mul(s, s, f);
                fmpz_mul(f, f, c);
                fmpz_mul(t, sum, xpow);
                fmpz_add(sum, t, s);
                fmpz_mod(sum, sum, mod);
            }
        }
    }

    /* Divide by factorial, XXX: Improve */
    {
        fmpq_t q;
        (q->num) = *sum;
        (q->den) = *f;
        padic_set_fmpq(res, q, ctx);
    }

    _fmpz_vec_clear(xpow, nxpow);
    fmpz_clear(c);
    fmpz_clear(s);
    fmpz_clear(t);
    fmpz_clear(f);
    fmpz_clear(sum);

    if (alloc)
        fmpz_clear(mod);
}

/*
    Sets \code{rop} to the exponential of \code{op}, reduced 
    in the given context.

    Assumptions:

        - The $p$-adic valuation of \code{op} is positive.
        - \code{op} is not zero modulo $p^N$.

    Supports aliasing.
 */
void _padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (ctx->N == 1)
    {
        padic_one(rop, ctx);
    }
    else if (ctx->N == 2)
    {
        padic_t one;

        _padic_init(one);
        _padic_one(one);
        padic_set(rop, op, ctx);
        padic_add(rop, rop, one, ctx);
        _padic_clear(one);
    }
    else
    {
        fmpz_t f, m, t, x;
        long i, k;

        fmpz_init(f);
        fmpz_init(m);
        fmpz_init(t);
        fmpz_init(x);

        padic_get_fmpz(x, op, ctx);

        i = _padic_exp_bound(padic_val(op), ctx->N, ctx->p) - 1;
        k = fmpz_fits_si(ctx->p) ? (i - 1) / (fmpz_get_si(ctx->p) - 1) : 0;

        fmpz_pow_ui(m, ctx->p, ctx->N + k);

        fmpz_one(padic_unit(rop));
        fmpz_one(f);

        for (i--; i >= 0; i--)
        {
            fmpz_mul_ui(f, f, i + 1);
            fmpz_mul(t, padic_unit(rop), x);
            fmpz_add(padic_unit(rop), f, t);

            fmpz_mod(f, f, m);
            fmpz_mod(padic_unit(rop), padic_unit(rop), m);
        }

        k = fmpz_remove(t, f, ctx->p);
        _padic_inv(f, t, ctx->p, ctx->N);
        fmpz_pow_ui(t, ctx->p, k);

        fmpz_divexact(padic_unit(rop), padic_unit(rop), t);
        fmpz_mul(padic_unit(rop), padic_unit(rop), f);
        padic_val(rop) = 0;
        padic_reduce(rop, ctx);

        fmpz_clear(f);
        fmpz_clear(m);
        fmpz_clear(t);
        fmpz_clear(x);
    }
}

int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    const long N = ctx->N;
    const long v = padic_val(op);

    if (fmpz_is_zero(padic_unit(op)))
    {
        padic_one(rop, ctx);
        return 1;
    }

    if (*(ctx->p) == 2L)
    {
        if (v <= 1)
        {
            return 0;
        }
        else if (v < N)
        {
            /*_padic_exp(rop, op, ctx);*/
            padic_exp_bsplit(rop, op, ctx);
            return 1;
        }
        else
        {
            padic_one(rop, ctx);
            return 1;
        }
    }
    else
    {
        if (v <= 0)
        {
            return 0;
        }
        else if (v < N)
        {
            /*_padic_exp(rop, op, ctx);*/
            padic_exp_bsplit(rop, op, ctx);
            return 1;
        }
        else
        {
            padic_one(rop, ctx);
            return 1;
        }
    }
}

