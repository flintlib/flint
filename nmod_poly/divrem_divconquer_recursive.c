/*
    Copyright (C) 2008, 2009, 2011 William Hart
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
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_divrem_divconquer_recursive(mp_ptr Q, mp_ptr BQ, mp_ptr W, mp_ptr V,
                          mp_srcptr A, mp_srcptr B, slong lenB, nmod_t mod)
{
    if (lenB <= NMOD_DIVREM_DIVCONQUER_CUTOFF)
    {
        mp_ptr t = V;
        mp_ptr w = t + 2*lenB - 1;
        
        flint_mpn_copyi(t + lenB - 1, A + lenB - 1, lenB);
        flint_mpn_zero(t, lenB - 1);
        
        _nmod_poly_divrem_basecase(Q, BQ, w, t, 2 * lenB - 1, B, lenB, mod);
        
        /* BQ = A - R */
        _nmod_vec_neg(BQ, BQ, lenB - 1, mod);
    }
    else
    {
        const slong n2 = lenB / 2;
        const slong n1 = lenB - n2;

        mp_ptr W1 = W;
        mp_ptr W2 = W + n2;

        mp_srcptr p1 = A + 2 * n2;
        mp_srcptr p2;
        mp_srcptr d1 = B + n2;
        mp_srcptr d2 = B;
        mp_srcptr d3 = B + n1;
        mp_srcptr d4 = B;

        mp_ptr q1   = Q + n2;
        mp_ptr q2   = Q;
        mp_ptr dq1  = BQ + n2;
        mp_ptr d1q1 = BQ + n2 - (n1 - 1);

        mp_ptr d2q1, d3q2, d4q2, t;

        /* 
           Set q1 to p1 div d1, a 2 n1 - 1 by n1 division so q1 ends up 
           being of length n1;  low(d1q1) = d1 q1 is of length n1 - 1
         */

        _nmod_poly_divrem_divconquer_recursive(q1, d1q1, W1, V, p1, d1, n1, mod);

        /* 
           Compute bottom n1 + n2 - 1 coeffs of d2q1 = d2 q1
         */

        d2q1 = W1;
        _nmod_poly_mullow(d2q1, q1, n1, d2, n2, n1 + n2 - 1, mod);

        /* 
           Compute dq1 = d1 q1 x^n2 + d2 q1, of length n1 + n2 - 1
           Split it into a segment of length n1 - 1 at dq1 and a piece
           of length n2 at BQ.
         */

        flint_mpn_copyi(dq1, d2q1, n1 - 1);
        if (n2 > n1 - 1)
            BQ[0] = d2q1[n1 - 1];

        _nmod_vec_add(d1q1, d1q1, d2q1 + n2, n1 - 1, mod);

        /*
           Compute t = A/x^n2 - dq1, which has length 2 n1 + n2 - 1, but we 
           are not interested in the top n1 coeffs as they will be zero, so 
           this has effective length n1 + n2 - 1

           For the following division, we want to set {p2, 2 n2 - 1} to the 
           top 2 n2 - 1 coeffs of this

           Since the bottom n2 - 1 coeffs of p2 are irrelevant for the 
           division, we in fact set {t, n2} to the relevant coeffs
         */

        t = W1;
        _nmod_vec_sub(t, A + n2 + (n1 - 1), BQ, n2, mod);
        p2 = t - (n2 - 1);

        /*
           Compute q2 = t div d3, a 2 n2 - 1 by n2 division, so q2 will have 
           length n2; let low(d3q2) = d3 q2, of length n2 - 1
         */

        d3q2 = BQ;
        _nmod_poly_divrem_divconquer_recursive(q2, d3q2, W2, V, p2, d3, n2, mod);

        /*
           Compute d4q2 = d4 q2, of length n1 + n2 - 1
         */

        d4q2 = W1;
        _nmod_poly_mullow(d4q2, d4, n1, q2, n2, n1 + n2 - 1, mod);

        /*
           Compute dq2 = d3q2 x^n1 + d4q2, of length n1 + n2 - 1
         */

        _nmod_vec_add(BQ + n1, BQ + n1, d3q2, n2 - 1, mod);
        flint_mpn_copyi(BQ, d4q2, n2);
        _nmod_vec_add(BQ + n2, BQ + n2, d4q2 + n2, n1 - 1, mod);

        /*
           Note Q = q1 x^n2 + q2, and BQ = dq1 x^n2 + dq2
         */
    }
}
