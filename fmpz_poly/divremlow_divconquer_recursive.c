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

    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#define FLINT_DIVREMLOW_DIVCONQUER_CUTOFF  16

void
_fmpz_poly_divremlow_divconquer_recursive(fmpz * Q, fmpz * QB, 
                                          const fmpz * A, const fmpz * B, long lenB)
{
    if (lenB <= FLINT_DIVREMLOW_DIVCONQUER_CUTOFF)
    {
        _fmpz_poly_divrem_basecase(Q, QB, A, 2 * lenB - 1, B, lenB);
        _fmpz_vec_sub(QB, A, QB, lenB - 1);
    }
    else
    {
        const long n2 = lenB / 2;
        const long n1 = lenB - n2;

        const fmpz * p1 = A + 2 * n2;
        const fmpz * p2;
        const fmpz * d1 = B + n2;
        const fmpz * d2 = B;

        fmpz * q1 = Q + n2;
        fmpz * q2 = Q;

        /*
           Think of the top lenB coefficients of QB as temporary space; this 
           code will not depend on {QB, lenB - 1} and W being adjacent
         */

        fmpz * W = QB + (lenB - 1);

        fmpz *d1q1, *d2q1, *t, *d3q2, *d4q2;

        /*
           Set q1 to p1 div d1, a 2 n1 - 1 by n1 division, so q1 has length 
           at most n1; {W, n1 - 1} is d1 * q1 is truncated to length n1 - 1
         */

        _fmpz_poly_divremlow_divconquer_recursive(q1, W, p1, d1, n1);

        /* 
           W is of length lenB, but we only care about the bottom n1 - 1 
           coeffs, which we push up by n2 + 1, to the very top; we do this 
           manually here instead of via _fmpz_vec_swap() because the source 
           and destination arrays overlap
         */

        d1q1 = W + (n2 + 1);

        {
            long i;
            for (i = 0; i < n1 - 1; i++)
                fmpz_swap(d1q1 + i, W + i);
        }

        /*
           Compute d2q1 = d2 * q1, of length at most lenB - 1; we'll need the 
           top n2 coeffs for t and the bottom n1 - 1 coeffs for QB
         */

        d2q1 = QB;

        _fmpz_poly_mul(d2q1, q1, n1, d2, n2);

        /*
           Compute {t - (n2 - 1), 2 n2 - 1} to be the top 2 n2 - 1 coeffs of 

               A / x^n2 - (d1q1 x^n2 + d2q1).

           Note that actually the bottom n2 - 1 coeffs may be arbitrary
         */

        t = W + n1;

        if (n1 == n2)
            fmpz_zero(t);
        _fmpz_vec_add(t, t, d2q1 + (n1 - 1), n2);
        _fmpz_vec_neg(t, t, n2);
        _fmpz_vec_add(t, t, A + (n1 + n2 - 1), n2);

        p2 = t - (n2 - 1);

        /*
           Move {QB, n1 - 1} into the bottom coefficients of W, so that 
           we can use {QB, 2 n2 - 1} as space in the next division
         */

        _fmpz_vec_swap(QB, W, n1 - 1);

        /*
           Compute q2 = t div {B + n1}, a 2 n2 - 1 by n2 division
         */

        d3q2 = QB;

        _fmpz_poly_divremlow_divconquer_recursive(q2, d3q2, p2, B + n1, n2);

        _fmpz_vec_swap(QB + n1, d3q2, n2 - 1);

        if (lenB & 1L)
            fmpz_zero(QB + n2);
        _fmpz_vec_add(QB + n2, QB + n2, W, n1 - 1);

        /*
           Compute {d4q2, lenB - 1} as {B, n1} * {q2, n2}; then move the 
           bottom n2 coeffs of this into {QB, n2}, and add the top n1 - 1 
           coeffs to {QB + n2, n1 - 1}
         */

        d4q2 = W;

        _fmpz_poly_mul(d4q2, B, n1, q2, n2);

        _fmpz_vec_swap(QB, d4q2, n2);
        _fmpz_vec_add(QB + n2, QB + n2, d4q2 + n2, n1 - 1);
    }
}

