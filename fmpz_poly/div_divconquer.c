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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_div_divconquer(fmpz * Q, const fmpz * A, long lenA, 
                                    const fmpz * B, long lenB)
{
    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 q1 - 1 by q1 division
         */

        const long q1 = lenA - lenB + 1;
        const long q2 = lenB - q1;

        fmpz * temp = _fmpz_vec_init(2 * q1 - 1);

        _fmpz_poly_div_divconquer_recursive(Q, temp, A + q2, B + q2, q1);

        _fmpz_vec_clear(temp, 2 * q1 - 1);
    }
    else if (lenA > 2 * lenB - 1)
    {
        /*
           Shift A right until it is of length 2 lenB - 1, call this p1
         */

        const long shift = lenA - 2 * lenB + 1;
        const fmpz * p1  = A + shift;

        fmpz * q2   = Q;
        fmpz * q1   = Q + shift;
        fmpz * dq1  = _fmpz_vec_init(lenA);
        fmpz * d1q1 = dq1 + shift;

        /* 
           Set q1 to p1 div B, a 2 lenB - 1 by lenB division, so q1 ends up 
           being of length lenB; also set d1q1 = d1 q1 of length 2 lenB - 1, 
           truncated to length lenB - 1
         */

        _fmpz_poly_divremlow_divconquer_recursive(q1, d1q1, p1, B, lenB);

        /* 
           We have dq1 = d1 q1 x^shift, of length lenA

           Compute R = A - dq1;  the first lenB coefficients represent 
           remainder terms (zero if division exact), leaving lenA - lenB 
           significant terms which we use in the division
         */

        _fmpz_vec_copy(dq1, A, shift);
        _fmpz_vec_sub(dq1 + shift, A + shift, dq1 + shift, lenB - 1);

        /*
           Compute q2 = trunc(R) div B;  it is a smaller division than the 
           original since len(trunc(R)) = lenA - lenB
         */

        _fmpz_poly_div_divconquer(q2, dq1, lenA - lenB, B, lenB);

        /*
           We have Q = q1 x^shift + q2;  Q has length lenB + shift, q2 has 
           length shift since it is an lenA - lenB by lenB division
         */

        _fmpz_vec_clear(dq1, lenA);
    }
    else  /* lenA = 2 lenB - 1 */
    {
        fmpz * temp = _fmpz_vec_init(lenA);

        _fmpz_poly_div_divconquer_recursive(Q, temp, A, B, lenB);

        _fmpz_vec_clear(temp, lenA);
    }
}

void
fmpz_poly_div_divconquer(fmpz_poly_t Q, 
                         const fmpz_poly_t A, const fmpz_poly_t B)
{
    long lenA, lenB;
    fmpz_poly_t t;
    fmpz * q;

    lenA = A->length;
    lenB = B->length;

    if (lenB == 0)
    {
        printf("Exception: division by zero in fmpz_poly_div_divconquer\n");
        abort();
    }

    if (lenA < lenB)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_poly_init2(t, lenA - lenB + 1);
        q = t->coeffs;
    }
    else
    {
        fmpz_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    _fmpz_poly_div_divconquer(q, A->coeffs, lenA, B->coeffs, lenB);

    if (Q == A || Q == B)
    {
        _fmpz_poly_set_length(t, lenA - lenB + 1);
        fmpz_poly_swap(t, Q);
        fmpz_poly_clear(t);
    }
    else
        _fmpz_poly_set_length(Q, lenA - lenB + 1);

    _fmpz_poly_normalise(Q);
}
