/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qadic.h"

void _qadic_pow(fmpz *rop, const fmpz *op, slong len, const fmpz_t e,
                   const fmpz *a, const slong *j, slong lena,
                   const fmpz_t p)
{
    const slong d = j[lena - 1];

    if (fmpz_is_zero(e))
    {
        fmpz_one(rop);
        _fmpz_vec_zero(rop + 1, 2 * d - 1 - 1);
    }
    else if (fmpz_is_one(e))
    {
        _fmpz_vec_set(rop, op, len);
        _fmpz_vec_zero(rop + len, 2 * d - 1 - len);
    }
    else
    {
        ulong bit;
        fmpz *v = _fmpz_vec_init(2 * d - 1);
        fmpz *R, *S, *T;

        _fmpz_vec_zero(rop, 2 * d - 1);

        /*
           Set bits to the bitmask with a 1 one place lower than the msb of e
         */

        bit = fmpz_bits(e) - 2;

        /*
           Trial run without any polynomial arithmetic to determine the parity
           of the number of swaps;  then set R and S accordingly
         */

        {
            unsigned int swaps = 0U;
            ulong bit2 = bit;
            if (fmpz_tstbit(e, bit2))
                swaps = ~swaps;
            while (bit2--)
                if (!fmpz_tstbit(e, bit2))
                    swaps = ~swaps;

            if (swaps == 0U)
            {
                R = rop;
                S = v;
            }
            else
            {
                R = v;
                S = rop;
            }
        }

        /*
           We unroll the first step of the loop, referring to {op, len}
         */

        _fmpz_poly_sqr(R, op, len);
        _fmpz_mod_poly_reduce(R, 2 * len - 1, a, j, lena, p);

        if (fmpz_tstbit(e, bit))
        {
            _fmpz_poly_mul(S, R, d, op, len);
            _fmpz_mod_poly_reduce(S, d + len - 1, a, j, lena, p);
            T = R;
            R = S;
            S = T;
        }

        while (bit--)
        {
            if (fmpz_tstbit(e, bit))
            {
                _fmpz_poly_sqr(S, R, d);
                _fmpz_mod_poly_reduce(S, 2 * d - 1, a, j, lena, p);
                _fmpz_poly_mul(R, S, d, op, len);
                _fmpz_mod_poly_reduce(R, d + len - 1, a, j, lena, p);
            }
            else
            {
                _fmpz_poly_sqr(S, R, d);
                _fmpz_mod_poly_reduce(S, 2 * d - 1, a, j, lena, p);
                T = R;
                R = S;
                S = T;
            }
        }

        _fmpz_vec_clear(v, 2 * d - 1);
    }
}

void qadic_pow(qadic_t x, const qadic_t y, const fmpz_t e, const qadic_ctx_t ctx)
{
    const slong N = qadic_prec(x);

    if (fmpz_sgn(e) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (qadic_pow).  e < 0.\n");
    }

    if (fmpz_is_zero(e))
    {
        qadic_one(x);
    }
    else if (qadic_is_zero(y))
    {
        qadic_zero(x);
    }
    else
    {
        fmpz_t val;  /* N - e * val(y) */

        fmpz_init_set(val, e);
        fmpz_mul_si(val, val, y->val);

        if (fmpz_cmp_si(val, N) >= 0)
        {
            qadic_zero(x);
        }
        else if (fmpz_is_one(e))
        {
            qadic_set(x, y, ctx);
        }
        else
        {
            const slong d = qadic_ctx_degree(ctx);
            fmpz *t;
            fmpz_t pow;
            int alloc;

            alloc = _padic_ctx_pow_ui(pow, N - fmpz_get_si(val), &ctx->pctx);

            if (x == y)
            {
                t = _fmpz_vec_init(2 * d - 1);
            }
            else
            {
                padic_poly_fit_length(x, 2 * d - 1);
                t = x->coeffs;
            }

            _qadic_pow(t, y->coeffs, y->length, e, ctx->a, ctx->j, ctx->len, pow);
            x->val = fmpz_get_si(val);

            if (x == y)
            {
                _fmpz_vec_clear(x->coeffs, x->alloc);
                x->coeffs = t;
                x->alloc  = 2 * d - 1;
                x->length = d;
            }
            else
            {
                _padic_poly_set_length(x, d);
            }
            _padic_poly_normalise(x);

            if (alloc)
                fmpz_clear(pow);
        }
        fmpz_clear(val);
    }
}

