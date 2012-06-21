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

    Copyright (C) 2011 Jan Tuitman
    Copyright (C) 2011, 2012 Sebastian Pancratz

******************************************************************************/

#include "padic.h"

void _padic_inv_precompute(padic_inv_t S, const fmpz_t p, long N)
{
#define n    (S->n)
#define pow  (S->pow)
#define u    (S->u)

    long *a, i;
    fmpz *t;

    n = FLINT_CLOG2(N) + 1;

    /* Compute sequence of exponents */
    a = flint_malloc(n * sizeof(long));
    for (a[i = 0] = N; a[i] > 1; i++)
        a[i + 1] = (a[i] + 1) / 2;

    pow = _fmpz_vec_init(2 * n + 2);
    u   = pow + n;
    t   = pow + 2 * n;

    /* Compute powers of p */
    {
        fmpz_one(t);
        fmpz_set(pow + i, p);
    }
    for (i--; i >= 1; i--)
    {
        if (a[i] & 1L)
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
        if (a[i] & 1L)
            fmpz_mul(pow + i, t, pow + (i + 1));
        else
            fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
    }

    flint_free(a);

#undef n
#undef pow
#undef u
}

void _padic_inv_clear(padic_inv_t S)
{
    _fmpz_vec_clear(S->pow, 2 * S->n + 2);
}

void _padic_inv_precomp(fmpz_t rop, const fmpz_t op, padic_inv_t S)
{
#define n    (S->n)
#define pow  (S->pow)
#define u    (S->u)

    long i;
    fmpz *t;

    t = pow + 2 * n;

    /* Compute reduced units */
    {
        fmpz_mod(u + 0, op, pow + 0);
    }
    for (i = 1; i < n; i++)
    {
        fmpz_mod(u + i, u + (i - 1), pow + i);
    }

    /* Run Newton iteration */
    i = n - 1;
    {
        fmpz_invmod(rop, u + i, pow + i);
    }
    for (i--; i >= 0; i--)
    {
        fmpz_mul(t, rop, rop);
        fmpz_mul(t + 1, u + i, t);
        fmpz_mul_2exp(rop, rop, 1);
        fmpz_sub(rop, rop, t + 1);
        fmpz_mod(rop, rop, pow + i);
    }

    fmpz_mod(rop, rop, pow + 0);

#undef n
#undef pow
#undef u
}

void _padic_inv(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
{
    if (N == 1)
    {
        fmpz_invmod(rop, op, p);
    }
    else
    {
        padic_inv_t S;

        _padic_inv_precompute(S, p, N);
        _padic_inv_precomp(rop, op, S);
        _padic_inv_clear(S);
    }
}

void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (_padic_is_zero(op))
    {
        printf("Exception (padic_inv).  Zero is not invertible.\n");
        abort();
    }

    /*
        If x = u p^v has negative valuation with N <= -v then the 
        exact inverse of x is zero when reduced modulo $p^N$
     */
    if (ctx->N + padic_val(op) <= 0)
    {
        padic_zero(rop);
        return;
    }

    _padic_inv(padic_unit(rop), 
               padic_unit(op), ctx->p, ctx->N + padic_val(op));

    padic_val(rop) = - padic_val(op);
}

