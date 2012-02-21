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

    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "padic.h"

extern long _padic_exp_bound(long v, long N, const fmpz_t p);

/*
    Computes the sum 
    \begin{equation*}
    (a-1)! x^{1-a} \sum_{i=a}^{b-1} \frac{x^i}{i!}.
    \end{equation*}
    in the rational $(T, Q)$.

    Assumes that $1 \leq a < b$.

    If $a + 1 = b$, sets $P = x$, $Q = a$, and $T = x$.
    If $a + 2 = b$, sets $P = x^2$, $Q = a (a + 1)$, $T = x (a + 1) + x^2$.
    In general, sets 
    \begin{align*}
    P & = x^{b-a}, \\
    Q & = \frac{(b-1)!}{(a-1)!}, \\
    T & = (b-1)! x^{1-a} \sum_{i=a}^{b-1} \frac{x^i}{i!}.
    \end{align*}
 */

static void
_padic_exp_bsplit_series(fmpz_t P, fmpz_t Q, fmpz_t T,
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

        _padic_exp_bsplit_series(P, Q, T, x, a, m);
        _padic_exp_bsplit_series(PR, QR, TR, x, m, b);

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
    Assumes that $x$ is such that $\exp(x)$ converges.

    Assumes that $v = \ord_p(x)$ with $v < N$, 
    which also forces $N$ to positive.

    The result $y$ might not be reduced modulo $p^N$.

    Supports aliasing between $x$ and $y$.
 */

static void
_padic_exp_bsplit(fmpz_t y, const fmpz_t x, long v, const fmpz_t p, long N)
{
    const long n = _padic_exp_bound(v, N, p);

    if (n == 1)
    {
        fmpz_one(y);
    }
    else
    {
        fmpz_t P, Q, T;

        fmpz_init(P);
        fmpz_init(Q);
        fmpz_init(T);

        _padic_exp_bsplit_series(P, Q, T, x, 1, n);

        fmpz_add(T, T, Q);  /* (T,Q) := (T,Q) + 1 */

        /* Note exp(x) is a unit so val(T) == val(Q). */
        if (fmpz_remove(T, T, p))
            fmpz_remove(Q, Q, p);

        _padic_inv(Q, Q, p, N);
        fmpz_mul(y, T, Q);

        fmpz_clear(P);
        fmpz_clear(Q);
        fmpz_clear(T);
    }
}

void 
_padic_exp_balanced_2(padic_t y, const padic_t x, const padic_ctx_t ctx)
{
    fmpz_t r, t, pv, pw, pN;
    long v;

    fmpz_init(r);
    fmpz_init(t);
    fmpz_init(pv);
    fmpz_init(pw);
    fmpz_init(pN);

    fmpz_pow_ui(t, ctx->p, padic_val(x) - 2);
    fmpz_mul(t, t, padic_unit(x));

    fmpz_set_ui(pv, 4);
    fmpz_set_ui(pw, 2);
    fmpz_pow_ui(pN, ctx->p, ctx->N);

    fmpz_one(padic_unit(y));
    padic_val(y) = 0L;

    for (v = 3; v < (2 * ctx->N); v *= 2)
    {
        fmpz_mul(pw, pw, pv);       /* pw = p^w, w = v - 1 */
        fmpz_mul(pv, pv, pv);       /* pv = p^v            */

        fmpz_fdiv_qr(t, r, t, pv);  /* r = p^w (t % p^v)   */
        fmpz_mul(r, r, pw);

        _padic_exp_bsplit(r, r, v - 1, ctx->p, ctx->N);
        fmpz_mul(padic_unit(y), padic_unit(y), r);
        fmpz_mod(padic_unit(y), padic_unit(y), pN);
    }

    fmpz_clear(r);
    fmpz_clear(t);
    fmpz_clear(pv);
    fmpz_clear(pw);
    fmpz_clear(pN);
}

_padic_exp_balanced_p(padic_t y, const padic_t x, const padic_ctx_t ctx)
{
    fmpz_t r, t, pv, pw, pN;
    long v;

    fmpz_init(r);
    fmpz_init(t);
    fmpz_init(pv);
    fmpz_init(pw);
    fmpz_init(pN);

    fmpz_pow_ui(t, ctx->p, padic_val(x) - 1);
    fmpz_mul(t, t, padic_unit(x));

    fmpz_set(pv, ctx->p);
    fmpz_one(pw);
    fmpz_pow_ui(pN, ctx->p, ctx->N);

    fmpz_one(padic_unit(y));
    padic_val(y) = 0L;

    for (v = 2; v < (2 * ctx->N); v *= 2)
    {
        fmpz_mul(pw, pw, pv);       /* pw = p^w, w = v - 1 */
        fmpz_mul(pv, pv, pv);       /* pv = p^v            */

        fmpz_fdiv_qr(t, r, t, pv);  /* r = p^w (t % p^v)   */
        fmpz_mul(r, r, pw);

        _padic_exp_bsplit(r, r, v - 1, ctx->p, ctx->N);
        fmpz_mul(padic_unit(y), padic_unit(y), r);
        fmpz_mod(padic_unit(y), padic_unit(y), pN);
    }

    fmpz_clear(r);
    fmpz_clear(t);
    fmpz_clear(pv);
    fmpz_clear(pw);
    fmpz_clear(pN);
}

/*
    Assumes that the exponential series converges at $x \neq 0$, 
    and that $\ord_p(x) < N$.

    Supports aliasing between $x$ and $y$.

    TODO:  Take advantage of additional factors of $p$ in $x$.
 */

void _padic_exp_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (*(ctx->p) == 2L)
        _padic_exp_balanced_2(rop, op, ctx);
    else
        _padic_exp_balanced_p(rop, op, ctx);
}

int padic_exp_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    const long N  = ctx->N;
    const long v  = padic_val(op);
    const fmpz *p = ctx->p;

    if (fmpz_is_zero(padic_unit(op)))
    {
        padic_one(rop, ctx);
        return 1;
    }

    if ((*p == 2L && v <= 1) || (v <= 0))
    {
        return 0;
    }
    else
    {
        if (v < N)
            _padic_exp_balanced(rop, op, ctx);
        else
            padic_one(rop, ctx);
        return 1;
    }
}
