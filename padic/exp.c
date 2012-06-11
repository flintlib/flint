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
    Returns an integer $i$ such that for all $j \geq i$ we have 
    $\ord_p(x^j / j!) \geq N$, where $\ord_p(x) = v$.

    When $p$ is a word-sized prime, 
    returns $\ceil{\frac{(p-1)N - 1}{(p-1)v - 1}}$.
    Otherwise, returns $\ceil{N/v}$.

    Assumes that $v < N$.  Moreover, $v$ has to be at least $2$ or $1$, 
    depending on whether $p$ is $2$ or odd.
 */
long _padic_exp_bound(long v, long N, const fmpz_t p)
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

/*
    Sets \code{rop} to the exponential of \code{op}, 
    not necessarily reduced modulo $p^N$.

    Assumes that the exponential series converges at $x \neq 0$, 
    where $x = p^v u$, and that $v = \ord_p(x) < N$.

    Supports aliasing.
 */
void _padic_exp_naive(fmpz_t rop, const fmpz_t u, long v, 
                                  const fmpz_t p, long N)
{
    const long n = _padic_exp_bound(v, N, p); 

    if (n == 1)
    {
        fmpz_one(rop);
    }
    else if (n == 2)
    {
        fmpz_t f;

        fmpz_init(f);
        fmpz_pow_ui(f, p, v);
        fmpz_mul(rop, f, u);
        fmpz_add_ui(rop, rop, 1);
        fmpz_clear(f);
    }
    else
    {
        fmpz_t f, m, t, x;
        long i, k;

        fmpz_init(f);
        fmpz_init(m);
        fmpz_init(t);
        fmpz_init(x);

        fmpz_pow_ui(f, p, v);
        fmpz_mul(x, f, u);

        i = n - 1;
        k = fmpz_fits_si(p) ? (i - 1) / (fmpz_get_si(p) - 1) : 0;

        fmpz_pow_ui(m, p, N + k);

        fmpz_one(rop);
        fmpz_one(f);

        for (i--; i >= 0; i--)
        {
            fmpz_mul_ui(f, f, i + 1);
            fmpz_mul(t, rop, x);
            fmpz_add(rop, f, t);

            fmpz_mod(f, f, m);
            fmpz_mod(rop, rop, m);
        }

        k = fmpz_remove(t, f, p);
        _padic_inv(f, t, p, N);
        fmpz_pow_ui(t, p, k);

        fmpz_divexact(rop, rop, t);
        fmpz_mul(rop, rop, f);

        fmpz_clear(f);
        fmpz_clear(m);
        fmpz_clear(t);
        fmpz_clear(x);
    }
}

void _padic_exp(fmpz_t rop, const fmpz_t u, long v, const fmpz_t p, long N)
{
    if (N < 1024)
        _padic_exp_rectangular(rop, u, v, p, N);
    else
        _padic_exp_balanced(rop, u, v, p, N);
}

int padic_exp(padic_t rop, const padic_t op, const padic_ctx_t ctx)
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
        {
            _padic_exp(padic_unit(rop), padic_unit(op), padic_val(op), p, N);
            padic_val(rop) = 0;
            _padic_reduce(rop, ctx);
        }
        else
        {
            padic_one(rop, ctx);
        }
        return 1;
    }
}

