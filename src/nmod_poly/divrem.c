/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 William Hart
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod.h"
#include "nmod_poly.h"
#include "gr_poly.h"

void
_nmod_poly_divrem(nn_ptr Q, nn_ptr R, nn_srcptr A, slong lenA,
                                  nn_srcptr B, slong lenB, nmod_t mod)
{
    gr_ctx_t ctx;
    _gr_ctx_init_nmod(ctx, &mod);

    /* the nmod method (_gr_nmod_poly_divrem) does the algorithm selection */
    if (_gr_poly_divrem(Q, R, A, lenA, B, lenB, ctx) != GR_SUCCESS)
        flint_throw(FLINT_IMPINV, "Exception (_nmod_poly_divrem). "
            "Leading coefficient of divisor is not invertible.\n");
}

void nmod_poly_divrem(nmod_poly_t Q, nmod_poly_t R,
                      const nmod_poly_t A, const nmod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    nmod_poly_t tQ, tR;
    nn_ptr q, r;

    if (lenB == 0)
    {
        if (nmod_poly_modulus(B) == 1)
        {
            nmod_poly_set(Q, A);
            nmod_poly_zero(R);
            return;
        } else
        {
            flint_throw(FLINT_DIVZERO, "Exception (nmod_poly_divrem). Division by zero.");
        }
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

    _nmod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, A->mod);

    if (Q == A || Q == B)
    {
        nmod_poly_swap(Q, tQ);
        nmod_poly_clear(tQ);
    }

    if (R == B)
    {
        nmod_poly_swap(R, tR);
        nmod_poly_clear(tR);
    }

    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}
