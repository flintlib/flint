/*
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
#include "ulong_extras.h"

void _nmod_poly_divrem_newton(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA, 
                              mp_srcptr B, slong lenB, nmod_t mod)
{
    const slong lenQ = lenA - lenB + 1;
    
    _nmod_poly_div_newton(Q, A, lenA, B, lenB, mod);

    if (lenB > 1)
    {
        if (lenQ >= lenB - 1)
            _nmod_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, mod);
        else
            _nmod_poly_mullow(R, B, lenB - 1, Q, lenQ, lenB - 1, mod);

        _nmod_vec_sub(R, A, R, lenB - 1, mod);
    }
}

void nmod_poly_divrem_newton(nmod_poly_t Q, nmod_poly_t R, 
                             const nmod_poly_t A, const nmod_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    mp_ptr q, r;

    if (lenB == 0)
    {
        if (nmod_poly_modulus(B) == 1)
        {
            nmod_poly_set(Q, A);
            nmod_poly_zero(R);
            return;
        } else
        {
            flint_printf("Exception (nmod_poly_divrem_newton). Division by zero.\n");
            flint_abort();
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
        q = _nmod_vec_init(lenA - lenB + 1);
    }
    else
    {
        nmod_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }
    if (R == A || R == B)
    {
        r = _nmod_vec_init(lenB - 1);
    }
    else
    {
        nmod_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem_newton(q, r, A->coeffs, lenA,
                                   B->coeffs, lenB, B->mod);

    if (Q == A || Q == B)
    {
        _nmod_vec_clear(Q->coeffs);
        Q->coeffs = q;
        Q->alloc  = lenA - lenB + 1;
    }
    if (R == A || R == B)
    {
        _nmod_vec_clear(R->coeffs);
        R->coeffs = r;
        R->alloc  = lenB - 1;
    }
    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _nmod_poly_normalise(R);
}
