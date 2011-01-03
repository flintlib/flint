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

******************************************************************************/

#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

void
_nmod_poly_divrem_divconquer(mp_ptr Q, mp_ptr R, mp_srcptr A, long lenA, 
                             mp_srcptr B, long lenB, nmod_t mod)
{
    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 n1 - 1 by n1 division
         */

        const long n1 = lenA - lenB + 1;
        const long n2 = lenB - n1;

        mp_srcptr p1 = A + n2;
        mp_srcptr d1 = B + n2;
        mp_srcptr d2 = B;

        mp_ptr V = nmod_vec_init((n1 - 1) + lenB - 1 + NMOD_DIVREM_DC_ITCH(n1, mod));
        mp_ptr W = V + NMOD_DIVREM_DC_ITCH(n1, mod);

        mp_ptr d1q1 = R + n2;
        mp_ptr d2q1 = W;

        _nmod_poly_divrem_divconquer_recursive(Q, d1q1, W, V, p1, d1, n1, mod);

        /*
           Compute d2q1 = Q d2, of length lenB - 1
         */

        if (n1 >= n2)
            _nmod_poly_mul(d2q1, Q, n1, d2, n2, mod);
        else
            _nmod_poly_mul(d2q1, d2, n2, Q, n1, mod);

        /*
           Compute BQ = d1q1 * x^n1 + d2q1, of length lenB - 1; 
           then compute R = A - BQ
         */

        mpn_copyi(R, d2q1, n2);
        _nmod_vec_add(R + n2, R + n2, d2q1 + n2, n1 - 1, mod);
        _nmod_vec_sub(R, A, R, lenB - 1, mod);

        nmod_vec_free(V);
    }
    else if (lenA > 2 * lenB - 1)
    {
        /*
           We shift A right until it is of length 2 lenB - 1, call this p1
         */

        const long shift = lenA - 2 * lenB + 1;
        mp_srcptr p1 = A + shift;

        mp_ptr V = nmod_vec_init(lenA + (2 * lenB - 1) + NMOD_DIVREM_DC_ITCH(lenB, mod));
        mp_ptr W = V + NMOD_DIVREM_DC_ITCH(lenB, mod);

        mp_ptr q1   = Q + shift;
        mp_ptr q2   = Q;
        mp_ptr dq1  = W;
        mp_ptr d1q1 = dq1 + shift;

        /* 
           Set q1 to p1 div B, a 2 lenB - 1 by lenB division, so q1 ends up 
           being of length lenB; set d1q1 = d1 * q1 of length 2 lenB - 1
         */

        _nmod_poly_divrem_divconquer_recursive(q1, d1q1, R, V, p1, B, lenB, mod);

        /* 
           We have dq1 = d1 * q1 * x^shift, of length lenA

           Compute R = A - dq1; the first lenB coeffs represent remainder 
           terms (zero if division is exact), leaving lenA - lenB significant 
           terms which we use in the division
         */

        mpn_copyi(dq1, A, shift);
        _nmod_vec_sub(dq1 + shift, A + shift, dq1 + shift, lenB - 1, mod);
        
        /*
           Compute q2 = trunc(R) div B; it is a smaller division than the 
           original since len(trunc(R)) = lenA - lenB
         */

        _nmod_poly_divrem_divconquer(q2, R, dq1, lenA - lenB, B, lenB, mod);

        /*
           We have Q = q1 * x^shift + q2; Q has length lenB + shift; 
           note q2 has length shift since the above division is 
           lenA - lenB by lenB

           We've also written the remainder in place
         */

        nmod_vec_free(V);
    }
    else  /* lenA = 2 * lenB - 1 */
    {
        mp_ptr V = nmod_vec_init(lenB - 1 + NMOD_DIVREM_DC_ITCH(lenB, mod));
        mp_ptr W = V + NMOD_DIVREM_DC_ITCH(lenB, mod);
 
        _nmod_poly_divrem_divconquer_recursive(Q, R, W, V, A, B, lenB, mod);
        _nmod_vec_sub(R, A, R, lenB - 1, mod);

        nmod_vec_free(V);
    }
}

void
nmod_poly_divrem_divconquer(nmod_poly_t Q, nmod_poly_t R,
                            const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tQ, tR;
    mp_ptr q, r;
    long Alen, Blen;

    Blen = B->length;

    if (Blen == 0)
    {
        printf("Exception: division by zero in nmod_poly_divrem_divconquer\n");
        abort();
    }

    Alen = A->length;

    if (Alen < Blen)
    {
        nmod_poly_zero(Q);
        nmod_poly_set(R, A);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2(tQ, A->mod.n, Alen - Blen + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, Alen - Blen + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_poly_init2(tR, A->mod.n, Blen - 1);
        r = tR->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, Blen - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem_divconquer(q, r, A->coeffs, Alen,
                                       B->coeffs, Blen, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(tQ, Q);
        nmod_poly_clear(tQ);
    }
    
    Q->length = Alen - Blen + 1;

    if (R == A || R == B)
    {
        nmod_poly_swap(tR, R);
        nmod_poly_clear(tR);
    }
    
    R->length = Blen - 1;

    _nmod_poly_normalise(R);
}
