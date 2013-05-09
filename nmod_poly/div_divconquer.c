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
    Copyright (C) 2012 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

static void
__nmod_poly_div_divconquer(mp_ptr Q, mp_srcptr A, len_t lenA, 
                             mp_srcptr B, len_t lenB, nmod_t mod)
{
    if (lenA < 2 * lenB - 1)
    {
        /*
           Convert unbalanced division into a 2 n1 - 1 by n1 division
         */

        const len_t n1 = lenA - lenB + 1;
        const len_t n2 = lenB - n1;

        mp_srcptr p1 = A + n2;
        mp_srcptr d1 = B + n2;
        
        mp_ptr V = _nmod_vec_init(n1 - 1 + NMOD_DIVREM_DC_ITCH(n1, mod));
        mp_ptr W = V + NMOD_DIVREM_DC_ITCH(n1, mod);

        _nmod_poly_div_divconquer_recursive(Q, W, V, p1, d1, n1, mod);

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


/* needed due to partial overlap */
static void
_nmod_vec_sub_dec(mp_ptr a, mp_srcptr b, mp_srcptr c, len_t n, nmod_t mod)
{
    len_t i;

    for (i = n - 1; i >= 0; i--)
        a[i] = n_submod(b[i], c[i], mod.n);
}

void
_nmod_poly_div_divconquer(mp_ptr Q, mp_srcptr A, len_t lenA,
                             mp_srcptr B, len_t lenB, nmod_t mod)
{
    if (lenA <= 2 * lenB - 1)
    {
        __nmod_poly_div_divconquer(Q, A, lenA, B, lenB, mod);
    }
    else  /* lenA > 2 * lenB - 1 */
    {
        mp_ptr S, T, V, R;
        len_t shift, next, n = 2 * lenB - 1;

        S = _nmod_vec_init(2 * n + (lenB - 1) + NMOD_DIVREM_DC_ITCH(lenB, mod));
        T = S + n;
        R = T + n;
        V = R + (lenB - 1);

        shift = lenA - n;
        _nmod_vec_set(S, A + shift, n);

        while (lenA >= n)
        {
            shift = lenA - n;
            _nmod_poly_divrem_divconquer_recursive(Q + shift, T, R, V, S, B, lenB, mod);
            next = FLINT_MIN(lenB, shift);
            _nmod_vec_sub_dec(S + next, S, T, lenB - 1, mod);
            _nmod_vec_set(S, A + shift - next, next);
            lenA -= lenB;
        }

        if (lenA >= lenB)
            __nmod_poly_div_divconquer(Q, S, lenA, B, lenB, mod);

        _nmod_vec_clear(S);
    }
}


void
nmod_poly_div_divconquer(nmod_poly_t Q,
                            const nmod_poly_t A, const nmod_poly_t B)
{
    nmod_poly_t tQ;
    mp_ptr q;
    len_t Alen, Blen;

    Blen = B->length;

    if (Blen == 0)
    {
        printf("Exception (nmod_poly_div_divconquer). Division by zero.\n");
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
