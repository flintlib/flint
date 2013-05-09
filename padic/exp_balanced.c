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
                        const fmpz_t x, len_t a, len_t b)
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
        const len_t m = (a + b) / 2;

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
_padic_exp_bsplit(fmpz_t y, const fmpz_t x, len_t v, const fmpz_t p, len_t N)
{
    const len_t n = _padic_exp_bound(v, N, p);

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

void _padic_exp_balanced_2(fmpz_t rop, const fmpz_t xu, len_t xv, len_t N)
{
    const fmpz_t p = {2L};

    fmpz_t r, t;
    len_t w;

    fmpz_init(r);
    fmpz_init(t);

    w = 1;

    fmpz_mul_2exp(t, xu, xv);
    fmpz_fdiv_r_2exp(t, t, N);

    fmpz_one(rop);

    while (!fmpz_is_zero(t))
    {
        fmpz_fdiv_r_2exp(r, t, 2*w);
        fmpz_sub(t, t, r);

        if (!fmpz_is_zero(r))
        {
            _padic_exp_bsplit(r, r, w, p, N);
            fmpz_mul(rop, rop, r);
            fmpz_fdiv_r_2exp(rop, rop, N);
        }

        w *= 2;
    }

    fmpz_clear(r);
    fmpz_clear(t);
}

void _padic_exp_balanced_p(fmpz_t rop, const fmpz_t xu, len_t xv, 
                                       const fmpz_t p, len_t N)
{
    fmpz_t r, t, pw, pN;
    len_t w;

    fmpz_init(r);
    fmpz_init(t);
    fmpz_init(pw);
    fmpz_init(pN);

    fmpz_set(pw, p);
    fmpz_pow_ui(pN, p, N);
    w = 1;

    fmpz_pow_ui(t, p, xv);
    fmpz_mul(t, t, xu);
    fmpz_mod(t, t, pN);

    fmpz_one(rop);

    while (!fmpz_is_zero(t))
    {
        fmpz_mul(pw, pw, pw);

        fmpz_fdiv_r(r, t, pw);
        fmpz_sub(t, t, r);

        if (!fmpz_is_zero(r))
        {
            _padic_exp_bsplit(r, r, w, p, N);
            fmpz_mul(rop, rop, r);
            fmpz_mod(rop, rop, pN);
        }

        w *= 2;
    }

    fmpz_clear(r);
    fmpz_clear(t);
    fmpz_clear(pw);
    fmpz_clear(pN);
}

/*
    Assumes that the exponential series converges at $x \neq 0$, 
    and that $\ord_p(x) < N$.

    Supports aliasing between $x$ and $y$.

    TODO:  Take advantage of additional factors of $p$ in $x$.
 */

void _padic_exp_balanced(fmpz_t rop, const fmpz_t u, len_t v, 
                                     const fmpz_t p, len_t N)
{
    if (fmpz_equal_ui(p, 2))
        _padic_exp_balanced_2(rop, u, v, N);
    else
        _padic_exp_balanced_p(rop, u, v, p, N);
}

int padic_exp_balanced(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    const len_t N  = padic_prec(rop);
    const len_t v  = padic_val(op);
    const fmpz *p = ctx->p;

    if (padic_is_zero(op))
    {
        padic_one(rop);
        return 1;
    }

    if ((fmpz_equal_ui(p, 2) && v <= 1) || (v <= 0))
    {
        return 0;
    }
    else
    {
        if (v < N)
        {
            _padic_exp_balanced(padic_unit(rop), padic_unit(op), padic_val(op), p, N);
            padic_val(rop) = 0;
        }
        else
        {
            padic_one(rop);
        }
        return 1;
    }
}
