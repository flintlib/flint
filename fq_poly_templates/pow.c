/*
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_pow) (TEMPLATE(T, struct) * rop,
                        const TEMPLATE(T, struct) * op, slong len, ulong e,
                        const TEMPLATE(T, ctx_t) ctx)
{
    ulong bit = ~((~UWORD(0)) >> 1);
    slong rlen;
    slong alloc = (slong) e * (len - 1) + 1;
    TEMPLATE(T, struct) * v = _TEMPLATE(T, vec_init) (alloc, ctx);
    TEMPLATE(T, struct) * R, *S, *T;

    /*
       Set bits to the bitmask with a 1 one place lower than the msb of e
     */

    while ((bit & e) == UWORD(0))
        bit >>= 1;

    bit >>= 1;

    /*
       Trial run without any polynomial arithmetic to determine the parity 
       of the number of swaps;  then set R and S accordingly
     */

    {
        unsigned int swaps = 0U;
        ulong bit2 = bit;
        if ((bit2 & e))
            swaps = ~swaps;
        while (bit2 >>= 1)
            if ((bit2 & e) == UWORD(0))
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
       We unroll the first step of the loop, referring to {poly, len}
     */

    _TEMPLATE(T, poly_sqr) (R, op, len, ctx);
    rlen = 2 * len - 1;
    if ((bit & e))
    {
        _TEMPLATE(T, poly_mul) (S, R, rlen, op, len, ctx);
        rlen += len - 1;
        T = R;
        R = S;
        S = T;
    }

    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            _TEMPLATE(T, poly_sqr) (S, R, rlen, ctx);
            rlen += rlen - 1;
            _TEMPLATE(T, poly_mul) (R, S, rlen, op, len, ctx);
            rlen += len - 1;
        }
        else
        {
            _TEMPLATE(T, poly_sqr) (S, R, rlen, ctx);
            rlen += rlen - 1;
            T = R;
            R = S;
            S = T;
        }
    }

    _TEMPLATE(T, vec_clear) (v, alloc, ctx);
}

void
TEMPLATE(T, poly_pow) (TEMPLATE(T, poly_t) rop, const TEMPLATE(T, poly_t) op,
                       ulong e, const TEMPLATE(T, ctx_t) ctx)
{
    const slong len = op->length;

    if ((len < 2) | (e < UWORD(3)))
    {
        if (e == UWORD(0))
            TEMPLATE(T, poly_one) (rop, ctx);
        else if (len == 0)
            TEMPLATE(T, poly_zero) (rop, ctx);
        else if (len == 1)
        {
            fmpz_t f;
            fmpz_init_set_ui(f, e);

            TEMPLATE(T, poly_fit_length) (rop, 1, ctx);
            TEMPLATE(T, pow) (rop->coeffs + 0, op->coeffs + 0, f, ctx);
            _TEMPLATE(T, poly_set_length) (rop, 1, ctx);

            fmpz_clear(f);
        }
        else if (e == UWORD(1))
            TEMPLATE(T, poly_set) (rop, op, ctx);
        else                    /* e == UWORD(2) */
            TEMPLATE(T, poly_sqr) (rop, op, ctx);
    }
    else
    {
        const slong rlen = (slong) e * (len - 1) + 1;

        if (rop != op)
        {
            TEMPLATE(T, poly_fit_length) (rop, rlen, ctx);
            _TEMPLATE(T, poly_pow) (rop->coeffs, op->coeffs, len, e, ctx);
            _TEMPLATE(T, poly_set_length) (rop, rlen, ctx);
        }
        else
        {
            TEMPLATE(T, poly_t) t;
            TEMPLATE(T, poly_init2) (t, rlen, ctx);
            _TEMPLATE(T, poly_pow) (t->coeffs, op->coeffs, len, e, ctx);
            _TEMPLATE(T, poly_set_length) (t, rlen, ctx);
            TEMPLATE(T, poly_swap) (rop, t, ctx);
            TEMPLATE(T, poly_clear) (t, ctx);
        }
    }
}


#endif
