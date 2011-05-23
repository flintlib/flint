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
    Assumes that \code{op} is a unit modulo $p^N$.

    In the current implementation supports aliasing, but this 
    might change.
 */
void _padic_inv_hensel(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
{
    if (N == 1)
    {
        fmpz_invmod(rop, op, p);
    }
    else
    {
        long *a, i, len;
        fmpz *W, *pow, *u;

        for (i = 1; (1L << i) < N; i++) ;

        /* Compute sequence of exponents */
        a = malloc((i + 1) * sizeof(long));
        for (a[i = 0] = N; a[i] > 1; i++)
            a[i + 1] = (a[i] + 1) / 2;
        len = i + 1;

        W   = _fmpz_vec_init(2 + 2 * len);
        pow = W + 2;
        u   = W + 2 + len;

        /* Compute powers of p */
        {
            fmpz_set_ui(W, 1);
            fmpz_set(pow + i, p);
        }
        for (i--; i >= 1; i--)
        {
            if (a[i] & 1L)
            {
                fmpz_mul(pow + i, W, pow + (i + 1));
                fmpz_mul(W, W, W);
            }
            else
            {
                fmpz_mul(W, W, pow + (i + 1));
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
            }
        }
        {
            if (a[i] & 1L)
                fmpz_mul(pow + i, W, pow + (i + 1));
            else
                fmpz_mul(pow + i, pow + (i + 1), pow + (i + 1));
        }

        /* Compute reduced units */
        {
            fmpz_mod(u, op, pow);
        }
        for (i = 1; i < len; i++)
        {
            fmpz_mod(u + i, u + (i - 1), pow + i);
        }

        /* Run Newton iteration */
        i = len - 1;
        {
            fmpz_invmod(rop, u + i, pow + i);
        }
        for (i--; i >= 0; i--)
        {
            fmpz_mul(W, u + i, rop);
            fmpz_sub_ui(W, W, 1);
            fmpz_mul(W + 1, W, rop);
            fmpz_sub(rop, rop, W + 1);
            fmpz_mod(rop, rop, pow + i);
        }

        free(a);
        _fmpz_vec_clear(W, 2 + 2 * len);
    }
}

void padic_inv_hensel(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (_padic_is_zero(op))
    {
        printf("Exception (padic_inv_hensel).  Zero is not invertible.\n");
        abort();
    }

    /*
        If x = u p^v has negative valuation with N <= -v 
        then there is no inverse of x defined modulo p^N.
     */
    if (ctx->N + padic_val(op) <= 0)
    {
        padic_zero(rop, ctx);
        return;
    }

    _padic_inv_hensel(padic_unit(rop), padic_unit(op), ctx->p, ctx->N + padic_val(op));

    padic_val(rop) = - padic_val(op);
}

