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

    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include "generics.h"

void
_elem_poly_pow_ui(elem_ptr res, elem_srcptr poly, long len, ulong e, const ring_t ring)
{
    ulong bit = ~((~0UL) >> 1);
    long rlen;
    long alloc = (long) e * (len - 1) + 1;
    elem_ptr v, R, S, T;

    v = _elem_vec_init(alloc, ring);

    while ((bit & e) == 0UL)
        bit >>= 1;
    bit >>= 1;

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

    _elem_poly_mul(R, poly, len, poly, len, ring);
    rlen = 2 * len - 1;

    if ((bit & e))
    {
        _elem_poly_mul(S, R, rlen, poly, len, ring);
        rlen += len - 1;
        T = R;
        R = S;
        S = T;
    }
    
    while ((bit >>= 1))
    {
        if ((bit & e))
        {
            /* should be sqr */
            _elem_poly_mul(S, R, rlen, R, rlen, ring);
            rlen += rlen - 1;
            _elem_poly_mul(R, S, rlen, poly, len, ring);
            rlen += len - 1;
        }
        else
        {
            /* should be sqr */
            _elem_poly_mul(S, R, rlen, R, rlen, ring);
            rlen += rlen - 1;
            T = R;
            R = S;
            S = T;
        }
    }

    _elem_vec_clear(v, alloc, ring);
}

void
elem_poly_pow_ui(elem_poly_t res, const elem_poly_t poly, ulong exp, const ring_t ring)
{
    long len = poly->length;

    if (exp == 0)
    {
        elem_one(res, ring);
    }
    else if (len == 0)
    {
        elem_zero(res, ring);
    }
    else if (exp == 1)
    {
        elem_poly_set(res, poly, ring);
    }
    else if (exp == 2)
    {
        /* should be sqr */
        elem_poly_mul(res, poly, poly, ring);
    }
    else if (res == poly)
    {
        elem_poly_t tmp;
        elem_init(tmp, ring);
        elem_poly_pow_ui(tmp, poly, exp, ring);
        elem_poly_swap(res, tmp);
        elem_clear(tmp, ring);
    }
    else
    {
        long rlen = (long) exp * (len - 1) + 1;
        elem_poly_fit_length(res, rlen, ring);
        _elem_poly_pow_ui(res->coeffs, poly->coeffs, len, exp, RING_PARENT(ring));
        elem_poly_set_length(res, rlen, ring);
        elem_poly_normalise(res, ring);
    }
}
