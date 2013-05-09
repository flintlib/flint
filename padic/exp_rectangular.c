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
    Computes the sum $1 + x + x^2 / 2$ reduced modulo $p^N$, 
    where $x = p^v u$.

    Supports aliasing between \code{rop} and $u$.
 */

static void _padic_exp_small(fmpz_t rop, const fmpz_t u, len_t v, len_t n, 
                                         const fmpz_t p, const fmpz_t pN)
{
    if (n == 1)  /* rop = 1 */
    {
        fmpz_one(rop);
    }
    else if (n == 2)  /* rop = 1 + x */
    {
        fmpz_t f;

        fmpz_init(f);
        fmpz_pow_ui(f, p, v);
        fmpz_mul(rop, f, u);
        fmpz_add_ui(rop, rop, 1);
        fmpz_mod(rop, rop, pN);
        fmpz_clear(f);
    }
    else  /* n == 3, rop = 1 + x + x^2 / 2 */
    {
        fmpz_t f;

        fmpz_init(f);
        fmpz_pow_ui(f, p, v);
        fmpz_mul(rop, f, u);
        fmpz_mul(f, rop, rop);
        if (fmpz_is_odd(f))
            fmpz_add(f, f, pN);
        fmpz_fdiv_q_2exp(f, f, 1);
        fmpz_add(rop, rop, f);
        fmpz_add_ui(rop, rop, 1);
        fmpz_clear(f);
    }
}

void _padic_exp_rectangular(fmpz_t rop, const fmpz_t u, len_t v, 
                                        const fmpz_t p, len_t N)
{
    const len_t n = _padic_exp_bound(v, N, p);

    fmpz_t pN;

    fmpz_init(pN);
    fmpz_pow_ui(pN, p, N);

    if (n <= 3)
    {
        _padic_exp_small(rop, u, v, n, p, pN);
    }
    else
    {
        const len_t k = fmpz_fits_si(p) ? 
                       (n - 1 - 1) / (fmpz_get_si(p) - 1) : 0;

        len_t i, npows, nsums;

        fmpz_t c, f, s, t, sum, pNk;
        fmpz *pows;

        fmpz_init(pNk);
        fmpz_pow_ui(pNk, p, N + k);

        npows = n_sqrt(n);
        nsums = (n + npows - 1) / npows;

        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(s);
        fmpz_init(t);
        fmpz_init(sum);

        /* Compute pows;  pows[i] = x^i. */
        pows = _fmpz_vec_init(npows + 1);
        fmpz_one(pows + 0);
        fmpz_pow_ui(f, p, v);
        fmpz_mul(pows + 1, f, u);
        for (i = 2; i <= npows; i++)
        {
            fmpz_mul(pows + i, pows + i - 1, pows + 1);
            fmpz_mod(pows + i, pows + i, pNk);
        }

        fmpz_zero(sum);
        fmpz_one(f);

        for (i = nsums - 1; i >= 0; i--)
        {
            len_t lo = i * npows;
            len_t hi = FLINT_MIN(n - 1, lo + npows - 1);

            fmpz_zero(s);
            fmpz_one(c);

            for ( ; hi >= lo; hi--)
            {
                fmpz_addmul(s, pows + hi - lo, c);
                if (hi != 0)
                    fmpz_mul_ui(c, c, hi);
            }

            fmpz_mul(t, pows + npows, sum);
            fmpz_mul(sum, s, f);
            fmpz_add(sum, sum, t);
            fmpz_mod(sum, sum, pNk);

            fmpz_mul(f, f, c);
        }

        /* Divide by factorial, TODO: Improve */

        /* Note exp(x) is a unit so val(sum) == val(f) */
        if (fmpz_remove(sum, sum, p))
            fmpz_remove(f, f, p);

        _padic_inv(f, f, p, N);
        fmpz_mul(rop, sum, f);

        _fmpz_vec_clear(pows, npows + 1);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(s);
        fmpz_clear(t);
        fmpz_clear(sum);
        fmpz_clear(pNk);
    }

    fmpz_mod(rop, rop, pN);
    fmpz_clear(pN);
}

int padic_exp_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
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
            _padic_exp(padic_unit(rop), padic_unit(op), padic_val(op), p, N);
            padic_val(rop) = 0;
        }
        else
        {
            padic_one(rop);
        }
        return 1;
    }
}

