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
#include "fmpz_mod_vec.h"
#include "fmpz_mod_poly.h"

#define FMPZ_MOD_DIV_DIVCONQUER_CUTOFF  16 /* must not exceed FMPZ_MOD_DIVREM_DIVCONQUER_CUTOFF */

void
_fmpz_mod_poly_div_divconquer_recursive(fmpz * Q, fmpz * W,
                                const fmpz * A, const fmpz * B, slong lenB,
                                              const fmpz_t invB, const fmpz_t p)
{
    if (lenB <= FMPZ_MOD_DIV_DIVCONQUER_CUTOFF)
    {
        _fmpz_mod_poly_div_basecase(Q, W, A, 2 * lenB - 1, B, lenB, invB, p);
    }
    else
    {
        const slong n2 = lenB / 2;
        const slong n1 = lenB - n2;

        fmpz * W1 = W;
        fmpz * W2 = W + lenB;

        const fmpz * p1 = A + 2 * n2;
        const fmpz * p2;
        const fmpz * d1 = B + n2;
        const fmpz * d2 = B;
        const fmpz * d3 = B + n1;
        
        fmpz * q1   = Q + n2;
        fmpz * q2   = Q;
        fmpz * d1q1 = W2;

        fmpz * d2q1, * t;

        /* 
           Set q1 to p1 div d1, a 2 n1 - 1 by n1 division so q1 ends up 
           being of length n1;  low(d1q1) = d1 q1 is of length n1 - 1
         */

        _fmpz_mod_poly_divrem_divconquer_recursive(q1, d1q1, W1, p1, d1, n1, invB, p);
        /* 
           Compute bottom n1 + n2 - 1 coeffs of d2q1 = d2 q1
         */

        d2q1 = W1;
        _fmpz_mod_poly_mullow(d2q1, q1, n1, d2, n2, p, n1 + n2 - 1);

        /* 
           Compute dq1 = d1 q1 x^n2 + d2 q1, of length n1 + n2 - 1
           Split it into a segment of length n1 - 1 at which is ignored 
           and a piece of length n2 at BQ.
         */

        if (n2 > n1 - 1)
            fmpz_set(W1 + 0, d2q1 + n1 - 1);

        _fmpz_mod_poly_add(W1 + n2 - (n1 - 1), d1q1, n1 - 1, d2q1 + n2, n1 - 1, p);
        
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
        _fmpz_mod_poly_sub(t, A + n2 + (n1 - 1), n2, t, n2, p);
        p2 = t - (n2 - 1);

        /*
           Compute q2 = t div d3, a 2 n2 - 1 by n2 division, so q2 will have 
           length n2; 
         */
        _fmpz_mod_poly_div_divconquer_recursive(q2, W2, p2, d3, n2, invB, p);

        /*
           Note Q = q1 x^n2 + q2
         */
    }
}
