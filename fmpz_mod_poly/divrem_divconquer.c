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
    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod_poly.h"

static void
__fmpz_mod_poly_divrem_divconquer(fmpz * Q, fmpz * R, 
    const fmpz * A, long lenA, const fmpz * B, long lenB, 
    const fmpz_t invB, const fmpz_t p)
{
    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 n1 - 1 by n1 division
         */

        const long n1 = lenA - lenB + 1;
        const long n2 = lenB - n1;

        const fmpz * p1 = A + n2;
        const fmpz * d1 = B + n2;
        const fmpz * d2 = B;

        fmpz * W = _fmpz_vec_init((2 * n1 - 1) + lenB - 1);

        fmpz * d1q1 = R + n2;
        fmpz * d2q1 = W + (2 * n1 - 1);

        _fmpz_mod_poly_divrem_divconquer_recursive(Q, d1q1, W, p1, d1, n1, 
                                                   invB, p);

        /*
           Compute d2q1 = Q d2, of length lenB - 1
         */

        if (n1 >= n2)
            _fmpz_mod_poly_mul(d2q1, Q, n1, d2, n2, p);
        else
            _fmpz_mod_poly_mul(d2q1, d2, n2, Q, n1, p);

        /*
           Compute BQ = d1q1 * x^n1 + d2q1, of length lenB - 1; 
           then compute R = A - BQ
         */

        _fmpz_vec_swap(R, d2q1, n2);
        _fmpz_mod_poly_add(R + n2, R + n2, n1 - 1, d2q1 + n2, n1 - 1, p);
        _fmpz_mod_poly_sub(R, A, lenA, R, lenA, p);

        _fmpz_vec_clear(W, (2 * n1 - 1) + lenB - 1);
    }
    else if (lenA > 2 * lenB - 1)
    {
        /*
           We shift A right until it is of length 2 lenB - 1, call this p1
         */

        const long shift = lenA - 2 * lenB + 1;
        const fmpz * p1 = A + shift;

        fmpz * q1   = Q + shift;
        fmpz * q2   = Q;
        fmpz * W    = R + lenA;
        fmpz * d1q1 = W + (2 * lenB - 1);

        /*
           XXX:  In this case, we expect R to be of length 
           lenA + 2 * (2 * lenB - 1) and A to be modifiable
         */

        /* 
           Set q1 to p1 div B, a 2 lenB - 1 by lenB division, so q1 ends up 
           being of length lenB; set d1q1 = d1 * q1 of length 2 lenB - 1
         */

        _fmpz_mod_poly_divrem_divconquer_recursive(q1, d1q1, W, p1, B, lenB, 
                                                   invB, p);

        /* 
           We have dq1 = d1 * q1 * x^shift, of length lenA

           Compute R = A - dq1; the first lenB coeffs represent remainder 
           terms (zero if division is exact), leaving lenA - lenB significant 
           terms which we use in the division
         */

        _fmpz_mod_poly_sub((fmpz *) A + shift, 
                           A + shift, lenB - 1, d1q1, lenB - 1, p);

        /*
           Compute q2 = trunc(R) div B; it is a smaller division than the 
           original since len(trunc(R)) = lenA - lenB
         */

        __fmpz_mod_poly_divrem_divconquer(q2, R, 
                                          A, lenA - lenB, B, lenB, invB, p);

        _fmpz_mod_poly_sub(R + lenA - lenB, 
                           A + lenA - lenB, lenB, d1q1 + lenB - 1, lenB, p);

        /*
           We have Q = q1 * x^shift + q2; Q has length lenB + shift; 
           note q2 has length shift since the above division is 
           lenA - lenB by lenB

           We've also written the remainder in place
         */
    }
    else  /* lenA = 2 * lenB - 1 */
    {
        fmpz * W = _fmpz_vec_init(lenA);

        _fmpz_mod_poly_divrem_divconquer_recursive(Q, R, W, 
                                                   A, B, lenB, invB, p);

        _fmpz_mod_poly_sub(R, A, lenB - 1, R, lenB - 1, p);

        _fmpz_vec_clear(W, lenA);
    }
}

void _fmpz_mod_poly_divrem_divconquer(fmpz *Q, fmpz *R, 
    const fmpz *A, long lenA, const fmpz *B, long lenB, 
    const fmpz_t invB, const fmpz_t p)
{
    if (lenA <= 2 * lenB - 1)
    {
        __fmpz_mod_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, invB, p);
    }
    else  /* lenA > 2 * lenB - 1 */
    {
        fmpz *S, *T;

        S = _fmpz_vec_init(2 * lenA + 2 * (2 * lenB - 1));
        T = S + lenA;
        _fmpz_vec_set(S, A, lenA);

        __fmpz_mod_poly_divrem_divconquer(Q, T, S, lenA, B, lenB, invB, p);

        _fmpz_vec_set(R, T, lenB - 1);
        _fmpz_vec_clear(S, 2 * lenA + 2 * (2 * lenB - 1));
    }
}

void
fmpz_mod_poly_divrem_divconquer(fmpz_mod_poly_t Q, fmpz_mod_poly_t R,
    const fmpz_mod_poly_t A, const fmpz_mod_poly_t B)
{
    const long lenA = A->length;
    const long lenB = B->length;
    const long lenQ = lenA - lenB + 1;

    fmpz *q, *r;
    fmpz_t invB;

    if (lenA < lenB)
    {
        fmpz_mod_poly_set(R, A);
        fmpz_mod_poly_zero(Q);
        return;
    }

    fmpz_init(invB);
    fmpz_invmod(invB, fmpz_mod_poly_lead(B), &(B->p));

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenQ);
    }
    else
    {
        fmpz_mod_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = _fmpz_vec_init(lenA);
    }
    else
    {
        fmpz_mod_poly_fit_length(R, lenA);
        r = R->coeffs;
    }

    _fmpz_mod_poly_divrem_divconquer(q, r, A->coeffs, lenA, 
                                           B->coeffs, lenB, invB, &(B->p));

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenQ;
        Q->length = lenQ;
    }
    else
    {
        _fmpz_mod_poly_set_length(Q, lenQ);
    }

    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenA;
        R->length = lenA;
    }
    _fmpz_mod_poly_set_length(R, lenB - 1);
    _fmpz_mod_poly_normalise(R);

    fmpz_clear(invB);
}
