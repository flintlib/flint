/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2010, 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T
	
#include "templates.h"
	
void
_TEMPLATE(T, poly_pow_trunc_binexp) (TEMPLATE(T, struct) * res,
                              const TEMPLATE(T, struct) * poly, ulong e,
                                     slong trunc, const TEMPLATE(T, ctx_t) ctx)
{
    ulong bit = ~((~UWORD(0)) >> 1);
    TEMPLATE(T, struct) * v = _TEMPLATE(T, vec_init) (trunc, ctx);
    TEMPLATE(T, struct) * R, * S, * T;

    /*
       Set bits to the bitmask with a 1 one place `lower than the msb of e
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
            R = res;
            S = v;
        }
        else
        {
            R = v;
            S = res;
        }
    }
    
    /*
       We unroll the first step of the loop, referring to {poly, len}
     */
    
    _TEMPLATE(T, poly_mullow) (R, poly, trunc, poly, trunc, trunc, ctx);
    if ((bit & e))
    {
        _TEMPLATE(T, poly_mullow) (S, R, trunc, poly, trunc, trunc, ctx);
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            _TEMPLATE(T, poly_mullow) (S, R, trunc, R, trunc, trunc, ctx);
            _TEMPLATE(T, poly_mullow) (R, S, trunc, poly, trunc, trunc, ctx);
        }
        else
        {
            _TEMPLATE(T, poly_mullow) (S, R, trunc, R, trunc, trunc, ctx);
            T = R;
            R = S;
            S = T;
        }
    }
    
    _TEMPLATE(T, vec_clear) (v, trunc, ctx);
}

void
TEMPLATE(T, poly_pow_trunc_binexp) (TEMPLATE(T, poly_t) res, 
                               const TEMPLATE(T, poly_t) poly, ulong e,
                                     slong trunc, const TEMPLATE(T, ctx_t) ctx)
{
    const slong len = poly->length;
    TEMPLATE(T, struct) * p;
    int pcopy = 0;

    if (len < 2 || e < UWORD(3) || trunc == 0)
    {
        if (len == 0 || trunc == 0)
            TEMPLATE(T, poly_zero) (res, ctx);
        else if (len == 1)
        {
            TEMPLATE(T, poly_fit_length) (res, 1, ctx);
            TEMPLATE(T, pow_ui) (res->coeffs + 0, poly->coeffs + 0, e, ctx);
            _TEMPLATE(T, poly_set_length) (res, 1, ctx);
            _TEMPLATE(T, poly_normalise) (res, ctx);
        }
        else if (e == UWORD(0))
        {
            TEMPLATE(T, t) c;
            TEMPLATE(T, init) (c, ctx);
            TEMPLATE(T, set_ui) (c, 1, ctx);
            TEMPLATE(T, poly_set_coeff) (res, 0, c, ctx);
            _TEMPLATE(T, poly_set_length) (res, 1, ctx);
            _TEMPLATE(T, poly_normalise) (res, ctx);
            TEMPLATE(T, clear) (c, ctx);
        }
        else if (e == UWORD(1))
        {
            TEMPLATE(T, poly_set) (res, poly, ctx);
            TEMPLATE(T, poly_truncate) (res, trunc, ctx);
        }
        else  /* e == UWORD(2) */
            TEMPLATE(T, poly_mullow) (res, poly, poly, trunc, ctx);

        return;
    }

    if (poly->length < trunc)
    {
        p = _TEMPLATE(T, vec_init) (trunc, ctx);
        _TEMPLATE(T, vec_set) (p, poly->coeffs, poly->length, ctx);
        _TEMPLATE(T, vec_zero) (p + poly->length, trunc - poly->length, ctx);
        pcopy = 1;
    } else
        p = poly->coeffs;

    if (res != poly || pcopy)
    {
        TEMPLATE(T, poly_fit_length) (res, trunc, ctx);
        _TEMPLATE(T, poly_pow_trunc_binexp) (res->coeffs, p, e, trunc, ctx);
    }
    else
    {
        TEMPLATE(T, poly_t) t;
        TEMPLATE(T, poly_init2) (t, trunc, ctx);
        _TEMPLATE(T, poly_pow_trunc_binexp) (t->coeffs, p, e, trunc, ctx);
        TEMPLATE(T, poly_swap) (res, t, ctx);
        TEMPLATE(T, poly_clear) (t, ctx);
    }

    if (pcopy)
        _TEMPLATE(T, vec_clear) (p, trunc, ctx);

    res->length = trunc;
    _TEMPLATE(T, poly_normalise) (res, ctx);
}

#endif
