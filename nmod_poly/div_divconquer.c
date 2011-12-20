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
_nmod_poly_div_divconquer(mp_ptr Q, mp_srcptr A, long lenA, 
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
        
        mp_ptr V = _nmod_vec_init(n1 - 1 + NMOD_DIVREM_DC_ITCH(n1, mod));
        mp_ptr W = V + NMOD_DIVREM_DC_ITCH(n1, mod);

        _nmod_poly_div_divconquer_recursive(Q, W, V, p1, d1, n1, mod);

        _nmod_vec_clear(V);
    }
    else if (lenA > 2 * lenB - 1)
    {
        /*
           We shift A right until it is of length 2 lenB - 1, call this p1
         */

        const long shift = lenA - 2 * lenB + 1;
        mp_srcptr p1 = A + shift;

        mp_ptr V = _nmod_vec_init(lenA + (3 * lenB - 2) + NMOD_DIVREM_DC_ITCH(lenB, mod));
        mp_ptr R = V + NMOD_DIVREM_DC_ITCH(lenB, mod);
        mp_ptr W = R + lenB - 1;

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

        _nmod_poly_div_divconquer(q2, dq1, lenA - lenB, B, lenB, mod);

        /*
           We have Q = q1 * x^shift + q2; Q has length lenB + shift; 
           note q2 has length shift since the above division is 
           lenA - lenB by lenB
         */

        _nmod_vec_clear(V);
    }
    else  /* lenA = 2 * lenB - 1 */
    {
        mp_ptr V = _nmod_vec_init(lenB - 1 + NMOD_DIVREM_DC_ITCH(lenB, mod));
        mp_ptr W = V + NMOD_DIVREM_DC_ITCH(lenB, mod);
 
        _nmod_poly_div_divconquer_recursive(Q, W, V, A, B, lenB, mod);
        
        _nmod_vec_clear(V);
    }
}

void
nmod_poly_div_divconquer(nmod_poly_t Q,
                            const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tQ;
    mp_ptr q;
    long Alen, Blen;

    Blen = B->length;

    if (Blen == 0)
    {
        printf("Exception: division by zero in nmod_poly_div_divconquer\n");
        abort();
    }

    Alen = A->length;

    if (Alen < Blen)
    {
        nmod_poly_zero(Q);
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

    _nmod_poly_div_divconquer(q, A->coeffs, Alen,
                                       B->coeffs, Blen, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(tQ, Q);
        nmod_poly_clear(tQ);
    }
    
    Q->length = Alen - Blen + 1;
}
