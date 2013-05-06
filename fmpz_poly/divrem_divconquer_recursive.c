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

#define FLINT_DIVREM_DIVCONQUER_CUTOFF  16

void
_fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ, fmpz * W, 
                                       const fmpz * A, const fmpz * B,
                                       long lenB)
{
    if (lenB <= FLINT_DIVREM_DIVCONQUER_CUTOFF)
    {
        _fmpz_vec_zero(BQ, lenB - 1);
        _fmpz_vec_set(BQ + (lenB - 1), A + (lenB - 1), lenB);

        _fmpz_poly_divrem_basecase(Q, BQ, BQ, 2 * lenB - 1, B, lenB);

        _fmpz_vec_neg(BQ, BQ, lenB - 1);
        _fmpz_vec_sub(BQ + (lenB - 1), A + (lenB - 1), BQ + (lenB - 1), lenB);
    }
    else
    {
        const long n2 = lenB / 2;
        const long n1 = lenB - n2;

        fmpz * W1 = W;
        fmpz * W2 = W + lenB;

        const fmpz * p1 = A + 2 * n2;
        const fmpz * p2;
        const fmpz * d1 = B + n2;
        const fmpz * d2 = B;
        const fmpz * d3 = B + n1;
        const fmpz * d4 = B;

        fmpz * q1   = Q + n2;
        fmpz * q2   = Q;
        fmpz * dq1  = BQ + n2;
        fmpz * d1q1 = BQ + 2 * n2;

        fmpz *d2q1, *d3q2, *d4q2, *t;

        /* 
           Set q1 to p1 div d1, a 2 n1 - 1 by n1 division so q1 ends up 
           being of length n1;  d1q1 = d1 q1 is of length 2 n1 - 1
         */

        _fmpz_poly_divrem_divconquer_recursive(q1, d1q1, W1, p1, d1, n1);

        /* 
           Compute d2q1 = d2 q1, of length lenB - 1
         */

        d2q1 = W1;
        _fmpz_poly_mul(d2q1, q1, n1, d2, n2);

        /* 
           Compute dq1 = d1 q1 x^n2 + d2 q1, of length 2 n1 + n2 - 1
         */

        _fmpz_vec_swap(dq1, d2q1, n2);
        _fmpz_vec_add(dq1 + n2, dq1 + n2, d2q1 + n2, n1 - 1);

        /*
           Compute t = A/x^n2 - dq1, which has length 2 n1 + n2 - 1, but we 
           are not interested in the top n1 coeffs as they will be zero, so 
           this has effective length n1 + n2 - 1

           For the following division, we want to set {p2, 2 n2 - 1} to the 
           top 2 n2 - 1 coeffs of this

           Since the bottom n2 - 1 coeffs of p2 are irrelevant for the 
           division, we in fact set {t, n2} to the relevant coeffs
         */

        t = BQ;
        _fmpz_vec_sub(t, A + n2 + (n1 - 1), dq1 + (n1 - 1), n2);
        p2 = t - (n2 - 1);

        /*
           Compute q2 = t div d3, a 2 n2 - 1 by n2 division, so q2 will have 
           length n2; let d3q2 = d3 q2, of length 2 n2 - 1
         */

        d3q2 = W1;
        _fmpz_poly_divrem_divconquer_recursive(q2, d3q2, W2, p2, d3, n2);

        /*
           Compute d4q2 = d4 q2, of length n1 + n2 - 1 = lenB - 1
         */

        d4q2 = W2;
        _fmpz_poly_mul(d4q2, d4, n1, q2, n2);

        /*
           Compute dq2 = d3q2 x^n1 + d4q2, of length n1 + 2 n2 - 1
         */

        _fmpz_vec_swap(BQ, d4q2, n2);
        _fmpz_vec_add(BQ + n2, BQ + n2, d4q2 + n2, n1 - 1);
        _fmpz_vec_add(BQ + n1, BQ + n1, d3q2, 2 * n2 - 1);

        /*
           Note Q = q1 x^n2 + q2, and BQ = dq1 x^n2 + dq2
         */
    }
}
