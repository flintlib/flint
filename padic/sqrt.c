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

#include <assert.h>
#include "padic.h"

/*
    Returns whether \code{op} has a square root modulo $p^N$ and if 
    so sets \code{rop} to such an element.

    Assumes that \code{op} is a unit modulo $p^N$.  Assumes $p$ is an 
    odd prime.

    In the current implementation, allows aliasing.
 */
static int _padic_sqrt(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
{
    int ans;

    assert(fmpz_is_odd(p));

    if (N == 1)
    {
        ans = fmpz_sqrtmod(rop, op, p);
        return ans;
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

        W      = _fmpz_vec_init(2 + 2 * len);
        pow    = W + 2;
        u      = W + (2 + len);

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
            ans = fmpz_sqrtmod(rop, u + i, p);
            if (ans)
                fmpz_invmod(rop, rop, p);
            else
                goto exit;
        }
        for (i--; i >= 0; i--)  /* z := z - z (a z^2 - 1) / 2 */
        {
            fmpz_mul(W, rop, rop);
            fmpz_mul(W + 1, u + i, W);
            fmpz_sub_ui(W + 1, W + 1, 1);

            if (fmpz_is_odd(W + 1))
                fmpz_add(W + 1, W + 1, pow + i);
            fmpz_fdiv_q_2exp(W + 1, W + 1, 1);

            fmpz_mul(W, W + 1, rop);
            fmpz_sub(rop, rop, W);
            fmpz_mod(rop, rop, pow + i);
        }

        /* Invert modulo p^N */
        fmpz_mul(rop, rop, u);
        fmpz_mod(rop, rop, pow);

      exit:

        free(a);
        _fmpz_vec_clear(W, 2 + 2 * len);

        return ans;
    }
}

/*
    Returns whether \code{op} has a square root modulo $2^N$ and if 
    so sets \code{rop} to such an element.

    Assumes that \code{op} is a unit modulo $2^N$.

    In the current implementation, allows aliasing.
 */
static int _padic_sqrt2(fmpz_t rop, const fmpz_t op, long N)
{
    assert(fmpz_is_odd(op));

    if (fmpz_fdiv_ui(op, 8) != 1)
        return 0;

    if (N <= 3)
    {
        fmpz_set_ui(rop, 1);
    }
    else
    {
        long *a, i, len;
        fmpz *W, *pow, *u;

        for (i = 1; (1L << i) < N; i++) ;

        /* Compute sequence of exponents */
        a = malloc((i + 2) * sizeof(long));
        for (a[i = 0] = N; a[i] > 3; i++)
            a[i + 1] = (a[i] + 3) / 2;
        len = i + 1;

        W      = _fmpz_vec_init(2 + 2 * len);
        pow    = W + 2;
        u      = W + (2 + len);

        /* Compute powers of p */
        {
            fmpz_set_ui(W, 1);
            fmpz_set_ui(pow + i, 8);
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
                fmpz_mul(pow + i, W, pow + (i + 1));
                fmpz_mul_ui(pow + i, pow + i, 2);
                fmpz_mul(W, W, W);
                fmpz_mul_ui(W, W, 2);
            }
        }
        {
            if (a[i] & 1L)
            {
                fmpz_mul(pow + i, W, pow + (i + 1));
            }
            else
            {
                fmpz_mul(pow + i, W, pow + (i + 1));
                fmpz_mul_ui(pow + i, pow + i, 2);
            }
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
        fmpz_set_ui(rop, 1);
        for (i--; i >= 0; i--)  /* z := z - z (a z^2 - 1) / 2 */
        {
            fmpz_mul(W, rop, rop);
            fmpz_mul(W + 1, u + i, W);
            fmpz_sub_ui(W + 1, W + 1, 1);
            fmpz_fdiv_q_2exp(W + 1, W + 1, 1);
            fmpz_mul(W, W + 1, rop);
            fmpz_sub(rop, rop, W);
            fmpz_mod(rop, rop, pow + i);
        }

        /* Invert modulo p^N */
        fmpz_mul(rop, rop, u);
        fmpz_mod(rop, rop, pow);

        free(a);
        _fmpz_vec_clear(W, 2 + 2 * len);

    }
    return 1;
}

int padic_sqrt(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (_padic_is_zero(op))
    {
        padic_zero(rop, ctx);
        return 1;
    }
    if (padic_val(op) & 1L)
    {
        return 0;
    }

    padic_val(rop) = padic_val(op) / 2;

    if (padic_val(rop) >= ctx->N)
    {
        padic_zero(rop, ctx);
        return 1;
    }

    if (fmpz_cmp_ui(ctx->p, 2))
    {
        return _padic_sqrt(padic_unit(rop), 
                           padic_unit(op), ctx->p, ctx->N - padic_val(rop));
    }
    else
    {
        return _padic_sqrt2(padic_unit(rop), 
                            padic_unit(op), ctx->N - padic_val(rop));
    }
}

