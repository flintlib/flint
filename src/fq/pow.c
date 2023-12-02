/*
    Copyright (C) 2011, 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fq.h"

void
_fq_pow(fmpz * rop, const fmpz * op, slong len, const fmpz_t e,
        const fq_ctx_t ctx)
{
    const slong d = fq_ctx_degree(ctx);

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
        _fq_reduce(R, 2 * len - 1, ctx);

        if (fmpz_tstbit(e, bit))
        {
            _fmpz_poly_mul(S, R, d, op, len);
            _fq_reduce(S, d + len - 1, ctx);
            T = R;
            R = S;
            S = T;
        }

        while (bit--)
        {
            if (fmpz_tstbit(e, bit))
            {
                _fmpz_poly_sqr(S, R, d);
                _fq_reduce(S, 2 * d - 1, ctx);
                _fmpz_poly_mul(R, S, d, op, len);
                _fq_reduce(R, d + len - 1, ctx);
            }
            else
            {
                _fmpz_poly_sqr(S, R, d);
                _fq_reduce(S, 2 * d - 1, ctx);
                T = R;
                R = S;
                S = T;
            }
        }

        _fmpz_vec_clear(v, 2 * d - 1);
    }
}

void
fq_pow(fq_t rop, const fq_t op, const fmpz_t e, const fq_ctx_t ctx)
{
    if (fmpz_sgn(e) < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fq_pow).  e < 0.\n");
    }

    if (fmpz_is_zero(e))
    {
        fq_one(rop, ctx);
    }
    else if (fq_is_zero(op, ctx))
    {
        fq_zero(rop, ctx);
    }
    else if (fmpz_is_one(e))
    {
        fq_set(rop, op, ctx);
    }
    else
    {
        const slong d = fq_ctx_degree(ctx);
        fmpz *t;

        if (rop == op)
        {
            t = _fmpz_vec_init(2 * d - 1);
        }
        else
        {
            fmpz_poly_fit_length(rop, 2 * d - 1);
            t = rop->coeffs;
        }

        if (fmpz_cmpabs(e, fq_ctx_prime(ctx)) < 0)
            _fq_pow(t, op->coeffs, op->length, e, ctx);
        else
        {
            fmpz_t order, e_mod;
            fmpz_init(e_mod);
	    fmpz_init(order);
	    fq_ctx_order(order, ctx);
	    fmpz_sub_ui(order, order, 1);
            fmpz_mod(e_mod, e, order);
	    _fq_pow(t, op->coeffs, op->length, e_mod, ctx);
	    fmpz_clear(order);
            fmpz_clear(e_mod);
        }

        if (rop == op)
        {
            _fmpz_vec_clear(rop->coeffs, rop->alloc);
            rop->coeffs = t;
            rop->alloc = 2 * d - 1;
            rop->length = d;
        }
        else
        {
            _fmpz_poly_set_length(rop, d);
        }
        _fmpz_poly_normalise(rop);
    }
}
