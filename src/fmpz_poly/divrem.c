/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

int _fmpz_poly_divrem(fmpz * Q, fmpz * R,
     const fmpz * A, slong lenA, const fmpz * B, slong lenB, int exact)
{
    int res;

    if (lenB < 6)
        res = _fmpz_poly_divrem_basecase(Q, R, A, lenA, B, lenB, exact);
    else
        res = _fmpz_poly_divrem_divconquer(Q, R, A, lenA, B, lenB, exact);

    return res;
}

void fmpz_poly_divrem(fmpz_poly_t Q, fmpz_poly_t R,
                      const fmpz_poly_t A, const fmpz_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    fmpz *q, *r;

    if (lenB == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_divrem). Division by zero.\n");
    }

    if (lenA < lenB)
    {
        fmpz_poly_set(R, A);
        fmpz_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        q = _fmpz_vec_init(lenA - lenB + 1);
    }
    else
    {
        fmpz_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        r = _fmpz_vec_init(lenA);
    }
    else
    {
        fmpz_poly_fit_length(R, lenA);
        r = R->coeffs;
    }

    _fmpz_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, 0);

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenA - lenB + 1;
        Q->length = lenA - lenB + 1;
    }
    else
    {
        _fmpz_poly_set_length(Q, lenA - lenB + 1);
    }

    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenA;
        R->length = lenA;
    }
    else
    {
        _fmpz_poly_set_length(R, lenA);
    }

    _fmpz_poly_normalise(Q);
    _fmpz_poly_normalise(R);
}
