/*
    Copyright (C) 2011 William Hart
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

/* TODO: tighten this code */

void
_arb_poly_div(arb_ptr Q,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong prec)
{
    slong lenQ, lenB2;
    arb_ptr Arev, Brev;

    lenQ = lenA - lenB + 1;

    Arev = _arb_vec_init(2 * lenQ);
    Brev = Arev + lenQ;

    _arb_poly_reverse(Arev, A + (lenA - lenQ), lenQ, lenQ);

    if (lenB >= lenQ)
    {
        _arb_poly_reverse(Brev, B + (lenB - lenQ), lenQ, lenQ);
        lenB2 = lenQ;
    }
    else
    {
        _arb_poly_reverse(Brev, B, lenB, lenB);
        lenB2 = lenB;
    }

    _arb_poly_div_series(Q, Arev, lenQ, Brev, lenB2, lenQ, prec);
    _arb_poly_reverse(Q, Q, lenQ, lenQ);

    _arb_vec_clear(Arev, 2 * lenQ);
}

void _arb_poly_divrem(arb_ptr Q, arb_ptr R,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong prec)
{
    const slong lenQ = lenA - lenB + 1;
    _arb_poly_div(Q, A, lenA, B, lenB, prec);

    if (lenB > 1)
    {
        if (lenQ >= lenB - 1)
            _arb_poly_mullow(R, Q, lenQ, B, lenB - 1, lenB - 1, prec);
        else
            _arb_poly_mullow(R, B, lenB - 1, Q, lenQ, lenB - 1, prec);
        _arb_vec_sub(R, A, R, lenB - 1, prec);
    }
}

void _arb_poly_rem(arb_ptr R,
    arb_srcptr A, slong lenA,
    arb_srcptr B, slong lenB, slong prec)
{
    const slong lenQ = lenA - lenB + 1;
    arb_ptr Q = _arb_vec_init(lenQ);
    _arb_poly_divrem(Q, R, A, lenA, B, lenB, prec);
    _arb_vec_clear(Q, lenQ);
}

int arb_poly_divrem(arb_poly_t Q, arb_poly_t R,
                             const arb_poly_t A, const arb_poly_t B, slong prec)
{
    const slong lenA = A->length, lenB = B->length;

    if (lenB == 0 || arb_contains_zero(B->coeffs + lenB - 1))
    {
        return 0;
    }

    if (lenA < lenB)
    {
        arb_poly_set(R, A);
        arb_poly_zero(Q);
        return 1;
    }

    if (Q == A || Q == B)
    {
        arb_poly_t T;
        arb_poly_init(T);
        arb_poly_divrem(T, R, A, B, prec);
        arb_poly_swap(Q, T);
        arb_poly_clear(T);
        return 1;
    }

    if (R == A || R == B)
    {
        arb_poly_t U;
        arb_poly_init(U);
        arb_poly_divrem(Q, U, A, B, prec);
        arb_poly_swap(R, U);
        arb_poly_clear(U);
        return 1;
    }

    arb_poly_fit_length(Q, lenA - lenB + 1);
    arb_poly_fit_length(R, lenB - 1);

    _arb_poly_divrem(Q->coeffs, R->coeffs, A->coeffs, lenA,
                                   B->coeffs, lenB, prec);

    _arb_poly_set_length(Q, lenA - lenB + 1);
    _arb_poly_set_length(R, lenB - 1);
    _arb_poly_normalise(R);

    return 1;
}
