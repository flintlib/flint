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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include "padic.h"

/*
    Returns an integer $i$ such that for all $j \geq i$ 
    $\ord_p(x^j / j!) \geq N$, where $\ord_p(x) = v$.
 */
static long bound(long v, long N, const fmpz_t p)
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
    Assumptions:

        - The $p$-adic valuation of \code{op} is positive.
        - \code{op} is not zero modulo $p^N$.
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
        padic_add(rop, op, one, ctx);
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

        i = bound(padic_val(op), ctx->N, ctx->p) - 1;
        k = fmpz_fits_si(ctx->p) ? (i - 1) / (fmpz_get_si(ctx->p) - 1) : 0;

        fmpz_pow_ui(m, ctx->p, ctx->N + k);

        fmpz_set_ui(padic_unit(rop), 1);
        fmpz_set_ui(f, 1);

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
    if (padic_is_zero(op, ctx))
    {
        padic_one(rop, ctx);
        return 1;
    }

    if (*(ctx->p) == 2L)
    {
        if (padic_val(op) > 1)
        {
            _padic_exp(rop, op, ctx);
            return 1;
        }
        else
        {
            return 0;
        }
    }
    else
    {
        if (padic_val(op) > 0)
        {
            _padic_exp(rop, op, ctx);
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

