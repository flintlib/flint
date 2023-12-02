/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"

int
_fmpz_poly_div(fmpz * Q, const fmpz * A, slong lenA,
                                         const fmpz * B, slong lenB, int exact)
{
    return _fmpz_poly_div_divconquer(Q, A, lenA, B, lenB, exact);
}

void
fmpz_poly_div(fmpz_poly_t Q, const fmpz_poly_t A, const fmpz_poly_t B)
{
    fmpz_poly_t tQ;
    fmpz *q;

    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_div). Division by zero.\n");
    }

    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        fmpz_poly_init2(tQ, A->length - B->length + 1);
        q = tQ->coeffs;
    }
    else
    {
        fmpz_poly_fit_length(Q, A->length - B->length + 1);
        q = Q->coeffs;
    }

    _fmpz_poly_div(q, A->coeffs, A->length, B->coeffs, B->length, 0);

    if (Q == A || Q == B)
    {
        _fmpz_poly_set_length(tQ, A->length - B->length + 1);
        fmpz_poly_swap(tQ, Q);
        fmpz_poly_clear(tQ);
    }
    else
        _fmpz_poly_set_length(Q, A->length - B->length + 1);

    _fmpz_poly_normalise(Q);
}
