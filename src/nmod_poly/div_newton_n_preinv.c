/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2013 Martin Lee

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_vec.h"
#include "nmod_poly.h"

void _nmod_poly_div_newton_n_preinv(nn_ptr Q, nn_srcptr A, slong lenA,
        nn_srcptr FLINT_UNUSED(B), slong lenB, nn_srcptr Binv,
        slong lenBinv, nmod_t mod)
{
    const slong lenQ = lenA - lenB + 1;
    nn_ptr Arev;
    TMP_INIT;

    TMP_START;
    Arev = TMP_ALLOC(lenQ * sizeof(ulong));
    _nmod_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);
    _nmod_poly_mullow(Q, Arev, lenQ, Binv, FLINT_MIN(lenQ,lenBinv), lenQ, mod);
    _nmod_poly_reverse(Q, Q, lenQ, lenQ);
    TMP_END;
}

void nmod_poly_div_newton_n_preinv(nmod_poly_t Q, const nmod_poly_t A,
                                    const nmod_poly_t B, const nmod_poly_t Binv)
{
    const slong lenA = A->length, lenB = B->length, lenQ = lenA - lenB + 1,
               lenBinv = Binv->length;

    nn_ptr q;

    if (lenB == 0)
    {
        if (nmod_poly_modulus(B) == 1)
        {
            nmod_poly_set(Q, A);
            return;
        } else
        {
            flint_throw(FLINT_DIVZERO, "Exception (nmod_poly_div_newton_n_preinv). Division by zero.\n");
        }
    }

    if (lenA < lenB)
    {
        nmod_poly_zero(Q);
        return;
    }

    if (lenA > 2 * lenB - 2)
    {
        flint_printf ("Exception (nmod_poly_div_newton_n_preinv).\n");
    }

    if (Q == A || Q == B || Q == Binv)
    {
        q = (nn_ptr) flint_malloc(lenQ * sizeof(ulong));
    }
    else
    {
        nmod_poly_fit_length(Q, lenQ);
        q = Q->coeffs;
    }

    _nmod_poly_div_newton_n_preinv(q, A->coeffs, lenA, B->coeffs, lenB,
                                                Binv->coeffs, lenBinv, B->mod);

    if (Q == A || Q == B || Q == Binv)
    {
        flint_free(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = lenQ;
    }
    Q->length = lenQ;
}
