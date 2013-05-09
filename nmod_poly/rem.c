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

void _nmod_poly_rem(mp_ptr R, mp_srcptr A, len_t lenA, 
                              mp_srcptr B, len_t lenB, nmod_t mod)
{
    if (lenA - lenB == 1)
    {
        _nmod_poly_rem_q1(R, A, lenA, B, lenB, mod);
    }
    else if (lenA < NMOD_DIVREM_DIVCONQUER_CUTOFF)
    {
        mp_ptr W = _nmod_vec_init(NMOD_DIVREM_BC_ITCH(lenA, lenB, mod));

        _nmod_poly_rem_basecase(R, W, A, lenA, B, lenB, mod);
        _nmod_vec_clear(W);
    }
    else
    {
        mp_ptr Q = _nmod_vec_init(lenA - lenB + 1);

        _nmod_poly_divrem(Q, R, A, lenA, B, lenB, mod);
        _nmod_vec_clear(Q);
    }
}

void nmod_poly_rem(nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B)
{
    const len_t lenA = A->length, lenB = B->length;
    nmod_poly_t tR;
    mp_ptr r;

    if (lenB == 0)
    {
        printf("Exception (nmod_poly_rem). Division by zero.\n");
        abort();
    }
    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        return;
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

    _nmod_poly_rem(r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (R == A || R == B)
    {
        nmod_poly_swap(R, tR);
        nmod_poly_clear(tR);
    }
        
    R->length = lenB - 1;
    _nmod_poly_normalise(R);
}
