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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

#include "mpn_extras.h"

void
_nmod_poly_divrem(mp_ptr Q, mp_ptr R, mp_srcptr A, long lenA, 
                                  mp_srcptr B, long lenB, nmod_t mod)
{
    if (lenA == lenB)
        _nmod_poly_divrem_q0(Q, R, A, B, lenB, mod);
    else if (lenA == lenB + 1)
        _nmod_poly_divrem_q1(Q, R, A, lenA, B, lenB, mod);
    else if (lenB < 15)
    {
        mp_ptr W = _nmod_vec_init(NMOD_DIVREM_BC_ITCH(lenA, lenB, mod));
        _nmod_poly_divrem_basecase(Q, R, W, A, lenA, B, lenB, mod);
        _nmod_vec_clear(W);
    }
    else if (lenB < 6000)
        _nmod_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, mod);
    else
        _nmod_poly_divrem_newton(Q, R, A, lenA, B, lenB, mod);
}

void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R,
                      const nmod_poly_t A, const nmod_poly_t B)
{
    const long lenA = A->length, lenB = B->length;
    nmod_poly_t tQ, tR;
    mp_ptr q, r;
    
    if (lenB == 0)
    {
        printf("Exception (nmod_poly_divrem). Division by zero.");
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
        nmod_poly_init2_preinv(tQ, A->mod.n, A->mod.ninv, lenA - lenB + 1);
        q = tQ->coeffs;
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        nmod_poly_init2_preinv(tR, B->mod.n, B->mod.ninv, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, tQ);
        nmod_poly_clear(tQ);
    }
    if (R == A || R == B)
    {
        nmod_poly_swap(R, tR);
        nmod_poly_clear(tR);
    }
        
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}

