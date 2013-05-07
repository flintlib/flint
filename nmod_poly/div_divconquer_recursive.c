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

    Copyright (C) 2008, 2009, 2011 William Hart
    Copyright (C) 2010 Sebastian Pancratz
   
******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_div_divconquer_recursive(mp_ptr Q, mp_ptr W, mp_ptr V,
                          mp_srcptr A, mp_srcptr B, long lenB, nmod_t mod)
{
    if (lenB <= NMOD_DIV_DIVCONQUER_CUTOFF)
    {
        _nmod_poly_div_basecase(Q, V, A, 2 * lenB - 1, B, lenB, mod);
    }
    else
    {
        const long n2 = lenB / 2;
        const long n1 = lenB - n2;

        mp_ptr W1 = W;
        mp_ptr W2 = W + n2;

        mp_srcptr p1 = A + 2 * n2;
        mp_srcptr p2;
        mp_srcptr d1 = B + n2;
        mp_srcptr d2 = B;
        mp_srcptr d3 = B + n1;
        
        mp_ptr q1   = Q + n2;
        mp_ptr q2   = Q;
        mp_ptr d1q1 = q2 + n2 - (n1 - 1);

        mp_ptr d2q1, t;

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
           Split it into a segment of length n1 - 1 at which is ignored 
           and a piece of length n2 at BQ.
         */

        if (n2 > n1 - 1)
            W1[0] = d2q1[n1 - 1];

        _nmod_vec_add(W1 + n2 - (n1 - 1), d1q1, d2q1 + n2, n1 - 1, mod);
        
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
        _nmod_vec_sub(t, A + n2 + (n1 - 1), t, n2, mod);
        p2 = t - (n2 - 1);

        /*
           Compute q2 = t div d3, a 2 n2 - 1 by n2 division, so q2 will have 
           length n2; 
         */

        _nmod_poly_div_divconquer_recursive(q2, W2, V, p2, d3, n2, mod);

        /*
           Note Q = q1 x^n2 + q2
         */
    }
}
