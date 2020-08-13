/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart

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

void _nmod_poly_rem(mp_ptr R, mp_srcptr A, slong lenA, 
                              mp_srcptr B, slong lenB, nmod_t mod)
{
    TMP_INIT;

    if (lenA - lenB == 1)
    {
        _nmod_poly_rem_q1(R, A, lenA, B, lenB, mod);
    }
    else if (lenA < NMOD_DIVREM_DIVCONQUER_CUTOFF)
    {
        mp_ptr W;
        
        TMP_START;
        W = TMP_ALLOC(NMOD_DIVREM_BC_ITCH(lenA, lenB, mod)*sizeof(mp_limb_t));

        _nmod_poly_rem_basecase(R, W, A, lenA, B, lenB, mod);
        TMP_END;
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
    const slong lenA = A->length, lenB = B->length;
    nmod_poly_t tR;
    mp_ptr r;

    if (lenB == 0)
    {
        flint_printf("Exception (nmod_poly_rem). Division by zero.\n");
        flint_abort();
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
