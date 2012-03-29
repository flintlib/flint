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
    Assumes that the exponential series converges at $x \neq 0$, 
    and that $\ord_p(x) < N$.
 */

void 
_padic_exp_rectangular(padic_t res, const padic_t x, const padic_ctx_t ctx)
{
    const long nterms = _padic_exp_bound(padic_val(x), ctx->N, ctx->p);

    if (nterms <= 3)
    {
        _padic_exp_naive(res, x, ctx);
    }
    else
    {
        const long k = fmpz_fits_si(ctx->p) ? 
                       (nterms - 1 - 1) / (fmpz_get_si(ctx->p) - 1) : 0;

        long i, npows, nsums;

        fmpz_t mod;
        int alloc;

        fmpz_t c, f, s, t, sum;
        fmpz *pows;

        alloc = _padic_ctx_pow_ui(mod, ctx->N + k, ctx);

        npows = n_sqrt(nterms);
        nsums = (nterms + npows - 1) / npows;

        fmpz_init(c);
        fmpz_init(f);
        fmpz_init(s);
        fmpz_init(t);
        fmpz_init(sum);

        /* Compute pows;  pows[i] = x^i. */
        pows = _fmpz_vec_init(npows + 1);
        fmpz_one(pows + 0);
        padic_get_fmpz(pows + 1, x, ctx);
        for (i = 2; i <= npows; i++)
        {
            fmpz_mul(pows + i, pows + i - 1, pows + 1);
            fmpz_mod(pows + i, pows + i, mod);
        }

        fmpz_zero(sum);
        fmpz_one(f);

        for (i = nsums - 1; i >= 0; i--)
        {
            long lo = i * npows;
            long hi = FLINT_MIN(nterms - 1, lo + npows - 1);

            fmpz_set(s, pows + hi - lo);
            fmpz_set_ui(c, hi);
            hi--;

            for ( ; hi > lo; hi--)
            {
                fmpz_addmul(s, pows + hi - lo, c);
                fmpz_mul_ui(c, c, hi);
            }

            if (hi == lo)
            {
                fmpz_add(s, s, c);
                if (hi != 0)
                    fmpz_mul_ui(c, c, hi);
            }

            fmpz_mul(t, pows + npows, sum);
            fmpz_mul(sum, s, f);
            fmpz_add(sum, sum, t);
            fmpz_mod(sum, sum, mod);

            fmpz_mul(f, f, c);
        }

        /* Divide by factorial, TODO: Improve */

        /* Note exp(x) is a unit so val(sum) == val(f). */
        if (_fmpz_remove(sum, ctx->p, ctx->pinv))
            _fmpz_remove(f,   ctx->p, ctx->pinv);

        _padic_inv(f, f, ctx->p, ctx->N);
        fmpz_mul(padic_unit(res), sum, f);
        padic_val(res) = 0;
        _padic_reduce(res, ctx);

        _fmpz_vec_clear(pows, npows + 1);
        fmpz_clear(c);
        fmpz_clear(f);
        fmpz_clear(s);
        fmpz_clear(t);
        fmpz_clear(sum);
        if (alloc)
            fmpz_clear(mod);
    }
}

int padic_exp_rectangular(padic_t rop, const padic_t op, const padic_ctx_t ctx)
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
            _padic_exp_rectangular(rop, op, ctx);
        else
            padic_one(rop, ctx);
        return 1;
    }
}

