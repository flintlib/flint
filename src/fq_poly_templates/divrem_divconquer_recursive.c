/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2013 Mike Hansen

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifdef T

#include "templates.h"

void
_TEMPLATE(T, poly_divrem_divconquer_recursive) (
    TEMPLATE(T, struct) * Q,
    TEMPLATE(T, struct) * BQ,
    TEMPLATE(T, struct) * W,
    const TEMPLATE(T, struct) * A,
    const TEMPLATE(T, struct) * B, slong lenB,
    const TEMPLATE(T, t) invB,
    const TEMPLATE(T, ctx_t) ctx)
{
    if (lenB <= TEMPLATE(CAP_T, POLY_DIVREM_DIVCONQUER_CUTOFF))
    {
        _TEMPLATE(T, vec_zero) (BQ, lenB - 1, ctx);
        _TEMPLATE(T, vec_set) (BQ + (lenB - 1), A + (lenB - 1), lenB, ctx);

        _TEMPLATE(T, poly_divrem_basecase) (Q, BQ, BQ, 2 * lenB - 1, B, lenB,
                                            invB, ctx);

        _TEMPLATE(T, poly_neg) (BQ, BQ, lenB - 1, ctx);
        _TEMPLATE(T, vec_set) (BQ + (lenB - 1), A + (lenB - 1), lenB, ctx);
    }
    else
    {
        const slong n2 = lenB / 2;
        const slong n1 = lenB - n2;

        TEMPLATE(T, struct) * W1 = W;
        TEMPLATE(T, struct) * W2 = W + lenB;

        const TEMPLATE(T, struct) * p1 = A + 2 * n2;
        const TEMPLATE(T, struct) * p2;
        const TEMPLATE(T, struct) * d1 = B + n2;
        const TEMPLATE(T, struct) * d2 = B;
        const TEMPLATE(T, struct) * d3 = B + n1;
        const TEMPLATE(T, struct) * d4 = B;

        TEMPLATE(T, struct) * q1 = Q + n2;
        TEMPLATE(T, struct) * q2 = Q;
        TEMPLATE(T, struct) * dq1 = BQ + n2;
        TEMPLATE(T, struct) * d1q1 = BQ + 2 * n2;

        TEMPLATE(T, struct) * d2q1, *d3q2, *d4q2, *t;

        /* 
           Set q1 to p1 div d1, a 2 n1 - 1 by n1 division so q1 ends up 
           being of length n1;  d1q1 = d1 q1 is of length 2 n1 - 1
         */

        _TEMPLATE(T, poly_divrem_divconquer_recursive) (q1, d1q1, W1,
                                                        p1, d1, n1, invB, ctx);

        /* 
           Compute d2q1 = d2 q1, of length lenB - 1
         */

        d2q1 = W1;
        _TEMPLATE(T, poly_mul) (d2q1, q1, n1, d2, n2, ctx);

        /* 
           Compute dq1 = d1 q1 x^n2 + d2 q1, of length 2 n1 + n2 - 1
         */

        _TEMPLATE(T, vec_swap) (dq1, d2q1, n2, ctx);
        _TEMPLATE(T, poly_add) (dq1 + n2, dq1 + n2, n1 - 1, d2q1 + n2, n1 - 1,
                                ctx);

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
        _TEMPLATE(T, poly_sub) (t, A + n2 + (n1 - 1), n2, dq1 + (n1 - 1), n2,
                                ctx);
        p2 = t - (n2 - 1);

        /*
           Compute q2 = t div d3, a 2 n2 - 1 by n2 division, so q2 will have 
           length n2; let d3q2 = d3 q2, of length 2 n2 - 1
         */

        d3q2 = W1;
        _TEMPLATE(T, poly_divrem_divconquer_recursive) (q2, d3q2, W2,
                                                        p2, d3, n2, invB, ctx);

        /*
           Compute d4q2 = d4 q2, of length n1 + n2 - 1 = lenB - 1
         */

        d4q2 = W2;
        _TEMPLATE(T, poly_mul) (d4q2, d4, n1, q2, n2, ctx);

        /*
           Compute dq2 = d3q2 x^n1 + d4q2, of length n1 + 2 n2 - 1
         */

        _TEMPLATE(T, vec_swap) (BQ, d4q2, n2, ctx);
        _TEMPLATE(T, poly_add) (BQ + n2, BQ + n2, n1 - 1, d4q2 + n2, n1 - 1,
                                ctx);
        _TEMPLATE(T, poly_add) (BQ + n1, BQ + n1, 2 * n2 - 1, d3q2, 2 * n2 - 1,
                                ctx);

        /*
           Note Q = q1 x^n2 + q2, and BQ = dq1 x^n2 + dq2
         */
    }
}


#endif
