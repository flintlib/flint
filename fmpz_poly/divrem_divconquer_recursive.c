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
_fmpz_poly_divrem_divconquer_recursive(fmpz * Q, fmpz * BQ,
                                       const fmpz * A, const fmpz * B,
                                       long B_len)
{
    long crossover = 16;
    long A_len = 2 * B_len - 1, n1, n2;

    const fmpz *d1, *d2, *d3, *d4, *p1, *p2;
    fmpz *q1, *q2, *dq1, *dq2, *d1q1, *d2q1, *d3q2, *d4q2, *t, *W1, *W2;

    if (B_len <= crossover)
    {
        /*
           Use the classical algorithm to compute the
           quotient and remainder, then use A - R to compute BQ
         */

        _fmpz_poly_divrem_basecase(Q, BQ, A, A_len, B, B_len);
        _fmpz_vec_sub(BQ, A, BQ, A_len);
        return;
    }

    n1 = (B_len + 1) / 2;
    n2 = B_len - n1;

    /*
       To avoid repeated allocations and de-allocations, 
       we allocate all the space we need up front.
     */
    
    W1 = _fmpz_vec_init((n1 + n2 - 1) + (n1 + 2 * n2 - 1));
    W2 = W1 + (n1 + n2 - 1);

    /* We let B = d1*x^n2 + d2 */

    d1 = B + n2;
    d2 = B;
    d3 = B + n1;
    d4 = B;

    /* 
       We let A = a1*x^(n1+2*n2-1) + a2*x^(n1+n2-1) + a3 
       where a1 is length n1, a2 is length n2 and a3 is length n1+n2-1 
       We set p1 = a1*x^(n1-1)+ other terms, so it has 
       length 2*n1-1 
     */

    p1 = A + 2 * n2;

    /* 
       Set q1 to p1 div d1 
       This is a 2*n1-1 by n1 division so 
       q1 ends up being length n1
       d1q1 = d1*q1 is length 2*n1-1
     */

    dq1 = BQ + n2;
    d1q1 = dq1 + n2;
    q1 = Q + n2;
    _fmpz_poly_divrem_divconquer_recursive(q1, d1q1, p1, d1, n1);

    /* 
       Compute d2q1 = d2*q1 
       which ends up being length n1+n2-1
     */

    d2q1 = W1;
    _fmpz_poly_mul(d2q1, q1, n1, d2, n2);

    /* 
       Compute dq1 = d1*q1*x^n2 + d2*q1
       dq1 is then of length 2*n1+n2-1
     */

    _fmpz_vec_copy(dq1, d2q1, n2);
    _fmpz_vec_add(dq1 + n2, dq1 + n2, d2q1 + n2, n1 - 1);

    /*
       Compute t = A/x^n2 - dq1
       which has length 2*n1+n2-1, but we are not interested 
       in the first n1 coefficients as they will be zero, 
       so it has effective length n1+n2-1
       We set p2 to the top 2*n2-1 coefficients of this
     */

    t = W1;
    _fmpz_vec_sub(t, A + n2, dq1, n1 + n2 - 1);
    p2 = t + (n1 > n2);

    /*
       Compute q2 = t div d3
       It is an 2*n2-1 by n2 division, so
       the length of q2 will be n2
       Also compute d3q2 of length 2*n2-1
     */

    dq2 = W2;
    d3q2 = dq2 + n1;
    q2 = Q;
    _fmpz_poly_divrem_divconquer_recursive(q2, d3q2, p2, d3, n2);

    /*
       Compute d4q2 = d4*q2 which is of length 
       n1+n2-1
     */

    d4q2 = W1;
    _fmpz_poly_mul(d4q2, d4, n1, q2, n2);

    /*
       Compute dq2 = d3q2*x^n1 + d4q2
       which is of length n1+2*n2-1
     */

    _fmpz_vec_copy(dq2, d4q2, n1);
    _fmpz_vec_add(dq2 + n1, dq2 + n1, d4q2 + n1, n2 - 1);

    /*
       Note Q = q1*x^n2 + q2
       so Q has length n1 + n2

       Write out BQ = dq1*x^n2 + dq2
       BQ has length 2*(n1+n2)-1
     */

    _fmpz_vec_copy(BQ, dq2, n2);
    _fmpz_vec_add(BQ + n2, BQ + n2, dq2 + n2, n1 + n2 - 1);

    _fmpz_vec_clear(W1, n1 + n2 - 1);
}
