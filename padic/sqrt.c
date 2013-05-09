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

/*
    Returns whether \code{op} has a square root modulo $p^N$ and if 
    so sets \code{rop} to such an element.

    Assumes that \code{op} is a unit modulo $p^N$.  Assumes $p$ is an 
    odd prime.

    In the current implementation, allows aliasing.
 */
static int _padic_sqrt_p(fmpz_t rop, const fmpz_t op, const fmpz_t p, len_t N)
{
    int ans;

    if (N == 1)
    {
        ans = fmpz_sqrtmod(rop, op, p);
        return ans;
    }
    else
    {
        len_t *e, i, n;
        fmpz *W, *pow, *u;

        e = _padic_lifts_exps(&n, N);

        W   = _fmpz_vec_init(2 + 2 * n);
        pow = W + 2;
        u   = W + (2 + n);

        _padic_lifts_pows(pow, e, n, p);

        /* Compute reduced units */
        {
            fmpz_mod(u, op, pow);
        }
        for (i = 1; i < n; i++)
        {
            fmpz_mod(u + i, u + (i - 1), pow + i);
        }

        /*
            Run Newton iteration for the inverse square root, 
            using the update formula 
                z := z - z (u z^2 - 1) / 2
            for all but the last step.  The last step is 
            replaced with 
                b := u z                  mod p^{N'}
                z := b + z (u - b^2) / 2  mod p^{N}.
         */
        i = n - 1;
        {
            ans = fmpz_sqrtmod(rop, u + i, p);
            if (!ans)
                goto exit;
            fmpz_invmod(rop, rop, p);
        }
        for (i--; i >= 1; i--)
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
        {
            fmpz_mul(W, u + 1, rop);
            fmpz_mul(W + 1, W, W);
            fmpz_sub(W + 1, u + 0, W + 1);
            if (fmpz_is_odd(W + 1))
                fmpz_add(W + 1, W + 1, pow + 0);
            fmpz_fdiv_q_2exp(W + 1, W + 1, 1);
            fmpz_mul(rop, rop, W + 1);
            fmpz_add(rop, W, rop);
            fmpz_mod(rop, rop, pow + 0);
        }

      exit:

        flint_free(e);
        _fmpz_vec_clear(W, 2 + 2 * n);

        return ans;
    }
}

/*
    Returns whether \code{op} has a square root modulo $2^N$ and if 
    so sets \code{rop} to such an element.

    Assumes that \code{op} is a unit modulo $2^N$.

    In the current implementation, allows aliasing.
 */
static int _padic_sqrt_2(fmpz_t rop, const fmpz_t op, len_t N)
{
    if (fmpz_fdiv_ui(op, 8) != 1)
        return 0;

    if (N <= 3)
    {
        fmpz_one(rop);
    }
    else
    {
        len_t *e, i, n;
        fmpz *W, *u;

        i = FLINT_CLOG2(N);

        /* Compute sequence of exponents */
        e = flint_malloc((i + 2) * sizeof(len_t));
        for (e[i = 0] = N; e[i] > 3; i++)
            e[i + 1] = (e[i] + 3) / 2;
        n = i + 1;

        W = _fmpz_vec_init(2 + n);
        u = W + 2;

        /* Compute reduced units */
        {
            fmpz_fdiv_r_2exp(u, op, e[0]);
        }
        for (i = 1; i < n; i++)
        {
            fmpz_fdiv_r_2exp(u + i, u + (i - 1), e[i]);
        }

        /* Run Newton iteration */
        fmpz_one(rop);
        for (i = n - 2; i >= 1; i--)  /* z := z - z (a z^2 - 1) / 2 */
        {
            fmpz_mul(W, rop, rop);
            fmpz_mul(W + 1, u + i, W);
            fmpz_sub_ui(W + 1, W + 1, 1);
            fmpz_fdiv_q_2exp(W + 1, W + 1, 1);
            fmpz_mul(W, W + 1, rop);
            fmpz_sub(rop, rop, W);
            fmpz_fdiv_r_2exp(rop, rop, e[i]);
        }
        {
            fmpz_mul(W, u + 1, rop);
            fmpz_mul(W + 1, W, W);
            fmpz_sub(W + 1, u + 0, W + 1);
            fmpz_fdiv_q_2exp(W + 1, W + 1, 1);
            fmpz_mul(rop, rop, W + 1);
            fmpz_add(rop, W, rop);
        }
        fmpz_fdiv_r_2exp(rop, rop, e[0]);

        flint_free(e);
        _fmpz_vec_clear(W, 2 + n);

    }
    return 1;
}

int _padic_sqrt(fmpz_t rop, const fmpz_t op, const fmpz_t p, len_t N)
{
    if (fmpz_equal_ui(p, 2))
    {
        return _padic_sqrt_2(rop, op, N);
    }
    else
    {
        return _padic_sqrt_p(rop, op, p, N);
    }
}

int padic_sqrt(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op))
    {
        padic_zero(rop);
        return 1;
    }
    if (padic_val(op) & 1L)
    {
        return 0;
    }

    padic_val(rop) = padic_val(op) / 2;

    /*
        In this case, if there is a square root it will be 
        zero modulo $p^N$.  We only have to establish whether 
        or not the element \code{op} is a square.
     */
    if (padic_val(rop) >= padic_prec(rop))
    {
        int ans;

        if (fmpz_equal_ui(ctx->p, 2))
        {
            ans = (fmpz_fdiv_ui(padic_unit(op), 8) == 1);
        }
        else
        {
            ans = fmpz_sqrtmod(padic_unit(rop), padic_unit(op), ctx->p);
        }
        padic_zero(rop);

        return ans;
    }

    return _padic_sqrt(padic_unit(rop), 
                       padic_unit(op), ctx->p, padic_prec(rop) - padic_val(rop));
}

