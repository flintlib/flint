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
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include "padic.h"

/*
    Assumes that op is non-zero and that op has valuation v s.t. 
    v > -N.

    Does not support aliasing.
 */
void _padic_inv_hensel(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    long N = ctx->N + op[1];
    const fmpz *p = ctx->p;

    /* Unit part */
    if (N == 1)
    {
        fmpz_invmod(rop, op, p);
    }
    else
    {
        long *a, i;
        long K = 1;  /* TODO:  Choice of cutoff */
        fmpz *W, *pow;

        for (i = 1; (1L << i) < N; i++) ;

        a = malloc((i + 1) * sizeof(long));
        a[i = 0] = N;
        while (N > 1)
            a[++i] = (N = (N + 1) / 2);

        W   = _fmpz_vec_init(3);
        pow = W + 2;

        /* Base case */
        {
            fmpz_pow_ui(pow, p, K);

            fmpz_invmod(rop, op, pow);
        }

        for (i--; i >= 0; i--)
        {
            N = a[i];

            fmpz_mul(W, pow, pow);
            if (N & 1L)
                fmpz_divexact(pow, W, p);
            else
                fmpz_swap(pow, W);

            fmpz_mul(W, op, rop);
            fmpz_sub_ui(W, W, 1);
            fmpz_mod(W, W, pow);
            fmpz_mul(W + 1, W, rop);
            fmpz_sub(rop, rop, W + 1);
            fmpz_mod(rop, rop, pow);
        }

        _fmpz_vec_clear(W, 3);
        free(a);
    }

    /* Valuation part */
    rop[1] = -op[1];
}

void padic_inv_hensel(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx))
    {
        printf("Exception (padic_inv).  Zero is not invertible.\n");
        abort();
    }

    /*
        If x = u p^v has negative valuation with N <= -v 
        then there is no inverse of x defined modulo p^N.
     */
    if (ctx->N + op[1] <= 0)
    {
        padic_zero(rop, ctx);
        return;
    }

    if (rop != op)
    {
        _padic_inv_hensel(rop, op, ctx);
    }
    else
    {
        padic_t t;

        padic_init(t, ctx);
        _padic_inv_hensel(t, op, ctx);
        padic_swap(rop, t, ctx);
        padic_clear(t, ctx);
    }
}

