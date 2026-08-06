/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_poly.h"

void _nmod_poly_rem(nn_ptr R, nn_srcptr A, slong lenA,
                              nn_srcptr B, slong lenB, nmod_t mod)
{
    /* Special case for GCD */
    if (lenA - lenB == 1)
    {
        ulong Q[2];
        ulong invB;

        invB = (B[lenB - 1] == 1) ? 1 : n_invmod(B[lenB - 1], mod.n);
        _nmod_poly_divrem_q1_preinv1(Q, R, A, lenA, B, lenB, invB, mod);
    }
    else if (lenB >= 2)
    {
        nn_ptr Q;
        TMP_INIT;

        TMP_START;
        Q = TMP_ALLOC((lenA - lenB + 1) * sizeof(ulong));
        _nmod_poly_divrem(Q, R, A, lenA, B, lenB, mod);
        TMP_END;
    }
}

void nmod_poly_rem(nmod_poly_t R, const nmod_poly_t A, const nmod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    nmod_poly_t tR;
    nn_ptr r;

    if (lenB == 0)
    {
        flint_throw(FLINT_DIVZERO, "Exception (nmod_poly_rem). Division by zero.\n");
    }
    if (lenA < lenB)
    {
        nmod_poly_set(R, A);
        return;
    }

    if (R == B)
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

    if (R == B)
    {
        nmod_poly_swap(R, tR);
        nmod_poly_clear(tR);
    }

    R->length = lenB - 1;
    _nmod_poly_normalise(R);
}
