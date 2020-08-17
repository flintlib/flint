/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_pow_binexp(fmpz * res, const fmpz * poly, slong len, ulong e)
{
    ulong bit = ~((~UWORD(0)) >> 1);
    slong rlen;
    slong alloc = (slong) e * (len - 1) + 1;
    fmpz *v = _fmpz_vec_init(alloc);
    fmpz *R, *S, *T;

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
    
    _fmpz_poly_sqr(R, poly, len);
    rlen = 2 * len - 1;
    if ((bit & e))
    {
        _fmpz_poly_mul(S, R, rlen, poly, len);
        rlen += len - 1;
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            _fmpz_poly_sqr(S, R, rlen);
            rlen += rlen - 1;
            _fmpz_poly_mul(R, S, rlen, poly, len);
            rlen += len - 1;
        }
        else
        {
            _fmpz_poly_sqr(S, R, rlen);
            rlen += rlen - 1;
            T = R;
            R = S;
            S = T;
        }
    }
    
    _fmpz_vec_clear(v, alloc);
}

void
fmpz_poly_pow_binexp(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)
{
    const slong len = poly->length;
    slong rlen;

    if ((len < 2) | (e < UWORD(3)))
    {
        if (e == UWORD(0))
            fmpz_poly_set_ui(res, 1);
        else if (len == 0)
            fmpz_poly_zero(res);
        else if (len == 1)
        {
            fmpz_poly_fit_length(res, 1);
            fmpz_pow_ui(res->coeffs, poly->coeffs, e);
            _fmpz_poly_set_length(res, 1);
        }
        else if (e == UWORD(1))
            fmpz_poly_set(res, poly);
        else  /* e == UWORD(2) */
            fmpz_poly_sqr(res, poly);
        return;
    }

    rlen = (slong) e * (len - 1) + 1;

    if (res != poly)
    {
        fmpz_poly_fit_length(res, rlen);
        _fmpz_poly_pow_binexp(res->coeffs, poly->coeffs, len, e);
        _fmpz_poly_set_length(res, rlen);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, rlen);
        _fmpz_poly_pow_binexp(t->coeffs, poly->coeffs, len, e);
        _fmpz_poly_set_length(t, rlen);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
