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
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

static 
void __nmod_poly_divrem_divconquer(mp_ptr Q, mp_ptr R, 
                                   mp_srcptr A, long lenA, 
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

        mp_ptr V = _nmod_vec_init((n1 - 1) + lenB - 1 + NMOD_DIVREM_DC_ITCH(n1, mod));
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

        flint_mpn_copyi(R, d2q1, n2);
        _nmod_vec_add(R + n2, R + n2, d2q1 + n2, n1 - 1, mod);
        _nmod_vec_sub(R, A, R, lenB - 1, mod);

        _nmod_vec_clear(V);
    }
    else  /* lenA = 2 * lenB - 1 */
    {
        mp_ptr V = _nmod_vec_init(lenB - 1 + NMOD_DIVREM_DC_ITCH(lenB, mod));
        mp_ptr W = V + NMOD_DIVREM_DC_ITCH(lenB, mod);
 
        _nmod_poly_divrem_divconquer_recursive(Q, R, W, V, A, B, lenB, mod);
        _nmod_vec_sub(R, A, R, lenB - 1, mod);

        _nmod_vec_clear(V);
    }
}

void _nmod_poly_divrem_divconquer(mp_ptr Q, mp_ptr R, 
                                  mp_srcptr A, long lenA, 
                                  mp_srcptr B, long lenB, nmod_t mod)
{
    if (lenA <= 2 * lenB - 1)
    {
        __nmod_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, mod);
    }
    else  /* lenA > 2 * lenB - 1 */
    {
        long shift, n = 2 * lenB - 1;
        mp_ptr S, QB, W, V, T;

        S = _nmod_vec_init(lenA + 2 * (lenB - 1) + n + NMOD_DIVREM_DC_ITCH(lenB, mod));
        QB = S + lenA;
        W = QB + (lenB - 1);
        T = W + (lenB - 1);
        V = T + n;

        _nmod_vec_set(S, A, lenA);

        while (lenA >= n)
        {
            shift = lenA - n;
            _nmod_poly_divrem_divconquer_recursive(Q + shift, QB, W, V, S + shift, B, lenB, mod);
            _nmod_vec_sub(S + shift, S + shift, QB, lenB - 1, mod);
            lenA -= lenB;
        }

        if (lenA >= lenB)
        {
            __nmod_poly_divrem_divconquer(Q, T, S, lenA, B, lenB, mod);
            _nmod_vec_set(S, T, lenA);
        }

        _nmod_vec_set(R, S, lenB - 1);
        _nmod_vec_clear(S);
    }
}

void nmod_poly_divrem_divconquer(nmod_poly_t Q, nmod_poly_t R,
                                 const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tQ, tR;
    mp_ptr q, r;
    long lenA, lenB;

    lenA = A->length;
    lenB = B->length;

    if (lenB == 0)
    {
        printf("Exception (nmod_poly_divrem_divconquer). Division by zero.\n");
        abort();
    }

    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        nmod_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        nmod_poly_init2(tQ, A->mod.n, lenA - lenB + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_poly_init2(tR, A->mod.n, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem_divconquer(q, r, A->coeffs, lenA,
                                       B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(tQ, Q);
        nmod_poly_clear(tQ);
    }
    if (R == A || R == B)
    {
        nmod_poly_swap(tR, R);
        nmod_poly_clear(tR);
    }
    
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}
