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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_pow_binexp(fmpz * res, const fmpz * poly, long len, ulong e)
{
    ulong bit = ~((~0UL) >> 1);
    long rlen;
    long alloc = (long) e * (len - 1) + 1;
    fmpz *v = _fmpz_vec_init(alloc);
    fmpz *R, *S, *T;

    /*
       Set bits to the bitmask with a 1 one place lower than the msb of e
     */
    
    while ((bit & e) == 0UL)
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
            if ((bit2 & e) == 0UL)
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
    
    _fmpz_poly_mul(R, poly, len, poly, len);
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
            _fmpz_poly_mul(S, R, rlen, R, rlen);
            rlen += rlen - 1;
            _fmpz_poly_mul(R, S, rlen, poly, len);
            rlen += len - 1;
        }
        else
        {
            _fmpz_poly_mul(S, R, rlen, R, rlen);
            rlen += rlen - 1;
            T = R;
            R = S;
            S = T;
        }
    }
    
    _fmpz_vec_clear(v, alloc);
}

void
_fmpz_poly_pow_binexp_r2l(fmpz * res, const fmpz * poly, long len, ulong e)
{
    /*
       We need a couple of copies of the input data to allow aliasing in the 
       multiplication methods, call these A (= res), B, and C.  Two of them  
       will be occupied by the current result and the current power of the   
       base, giving rise to six positions pos for (A, B, C):                 
       0 - (res, base, empty)                                              
       1 - (res, empty, base)                                              
       2 - (base, res, empty)                                              
       3 - (base, empty, res)                                              
       4 - (empty, res, base)                                              
       5 - (empty, base, res)                                              
       We keep track of the length of the power of the base and the current  
       result in lenB and lenR.                                              
     */
    const long alloc = (long) e * (len - 1) + 1;
    fmpz *A, *B, *C;
    long lenB, lenR;
    int pos;

    A = res;
    B = _fmpz_vec_init(2 * alloc);
    C = B + alloc;
    _fmpz_vec_copy(B, poly, len);
    fmpz_set_ui(A, 1UL);
    lenB = len;
    lenR = 1;
    pos = 0;

    while (1)
    {
        if (e & 1UL)
        {
            switch (pos)
            {
                case 0:
                    _fmpz_poly_mul(C, B, lenB, A, lenR);
                    pos = 5;
                    break;
                case 1:
                    _fmpz_poly_mul(B, C, lenB, A, lenR);
                    pos = 4;
                    break;
                case 2:
                    _fmpz_poly_mul(C, A, lenB, B, lenR);
                    pos = 3;
                    break;
                case 3:
                    _fmpz_poly_mul(B, A, lenB, C, lenR);
                    pos = 2;
                    break;
                case 4:
                    _fmpz_poly_mul(A, C, lenB, B, lenR);
                    pos = 1;
                    break;
                case 5:
                    _fmpz_poly_mul(A, B, lenB, C, lenR);
                    pos = 0;
                    break;
            }
            lenR += lenB - 1;
        }

        e >>= 1;
        if (e == 0UL)
            break;

        switch (pos)
        {
            case 0:
                _fmpz_poly_mul(C, B, lenB, B, lenB);
                pos = 1;
                break;
            case 1:
                _fmpz_poly_mul(B, C, lenB, C, lenB);
                pos = 0;
                break;
            case 2:
                _fmpz_poly_mul(C, A, lenB, A, lenB);
                pos = 4;
                break;
            case 3:
                _fmpz_poly_mul(B, A, lenB, A, lenB);
                pos = 5;
                break;
            case 4:
                _fmpz_poly_mul(A, C, lenB, C, lenB);
                pos = 2;
                break;
            case 5:
                _fmpz_poly_mul(A, B, lenB, B, lenB);
                pos = 3;
                break;
        }
        lenB += lenB - 1;
    }

    if (pos == 2 || pos == 4)
        _fmpz_vec_swap(A, B, lenR);
    if (pos == 3 || pos == 5)
        _fmpz_vec_swap(A, C, lenR);

    _fmpz_vec_clear(B, alloc);
}

void
fmpz_poly_pow_binexp(fmpz_poly_t res, const fmpz_poly_t poly, ulong e)
{
    const long len = poly->length;
    long rlen;

    if ((len < 2) | (e < 3UL))
    {
        if (e == 0UL)
            fmpz_poly_set_ui(res, 1);
        else if (len == 0)
            fmpz_poly_zero(res);
        else if (len == 1)
        {
            fmpz_poly_fit_length(res, 1);
            fmpz_pow_ui(res->coeffs, poly->coeffs, e);
            _fmpz_poly_set_length(res, 1);
        }
        else if (e == 1UL)
            fmpz_poly_set(res, poly);
        else  /* e == 2UL */
            fmpz_poly_mul(res, poly, poly);
        return;
    }

    rlen = (long) e * (len - 1) + 1;

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
