/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_pow_trunc(fmpz * res, const fmpz * poly, ulong e, slong n)
{
    ulong bit = ~((~UWORD(0)) >> 1);
    fmpz *v = _fmpz_vec_init(n);
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
       We unroll the first step of the loop, referring to {poly, n}
     */
    
    _fmpz_poly_sqrlow(R, poly, n, n);
    if ((bit & e))
    {
        _fmpz_poly_mullow(S, R, n, poly, n, n);
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            _fmpz_poly_sqrlow(S, R, n, n);
            _fmpz_poly_mullow(R, S, n, poly, n, n);
        }
        else
        {
            _fmpz_poly_sqrlow(S, R, n, n);
            T = R;
            R = S;
            S = T;
        }
    }
    
    _fmpz_vec_clear(v, n);
}

void
fmpz_poly_pow_trunc(fmpz_poly_t res, const fmpz_poly_t poly, ulong e, slong n)
{
    fmpz * copy;
    int clear;
    slong i, len;

    if (n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }
    if (e == 0)
    {
        fmpz_poly_set_ui(res, 1);
        return;
    }

    /* Set len to the length of poly mod x^n */
    len = FLINT_MIN(n, poly->length);
    for (--len; (len >= 0) && !poly->coeffs[len]; --len) ;
    ++len;

    if ((len < 2) | (e < 3))
    {
        if (len == 0)
            fmpz_poly_zero(res);
        else if (len == 1)
        {
            fmpz_poly_fit_length(res, 1);
            fmpz_pow_ui(res->coeffs, poly->coeffs, e);
            _fmpz_poly_set_length(res, 1);
        }
        else if (e == 1)
        {
            if (res != poly)
            {
                fmpz_poly_fit_length(res, len);
                _fmpz_vec_set(res->coeffs, poly->coeffs, len);
                _fmpz_poly_set_length(res, len);
            }
            else
                fmpz_poly_truncate(res, len);
        }
        else  /* e == 2 */
            fmpz_poly_sqrlow(res, poly, n);
        return;
    }

    if (poly->length >= n)
    {
        copy = poly->coeffs;
        clear = 0;
    }
    else
    {
        copy = (fmpz *) flint_malloc(n * sizeof(fmpz));
        for (i = 0; i < poly->length; i++)
            copy[i] = poly->coeffs[i];
        flint_mpn_zero((mp_ptr) copy + poly->length, n - poly->length);
        clear = 1;
    }

    if (res != poly)
    {
        fmpz_poly_fit_length(res, n);
        _fmpz_poly_pow_trunc(res->coeffs, copy, e, n);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_pow_trunc(t->coeffs, copy, e, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);

    if (clear)
        flint_free(copy);
}

