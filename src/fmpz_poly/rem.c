/*
    Copyright (C) 2010, 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_rem(fmpz * R, const fmpz * A, slong lenA, const fmpz * B, slong lenB)
{
    if (lenA < 15)
    {
        _fmpz_poly_rem_basecase(R, A, lenA, B, lenB);
    }
    else
    {
        fmpz *Q = _fmpz_vec_init(lenA - lenB + 1);

        _fmpz_poly_divrem(Q, R, A, lenA, B, lenB, 0);

        _fmpz_vec_clear(Q, lenA - lenB + 1);
    }
}

void fmpz_poly_rem(fmpz_poly_t R, const fmpz_poly_t A, const fmpz_poly_t B)
{
    const slong lenA = A->length, lenB = B->length;
    fmpz *r;

    if (lenB == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_rem). Division by zero.\n");
    }

    if (lenA < lenB)
    {
        fmpz_poly_set(R, A);
        return;
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

    _fmpz_poly_rem(r, A->coeffs, lenA, B->coeffs, lenB);

    if (R == A || R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = lenB - 1;
    }
    _fmpz_poly_set_length(R, lenA);
    _fmpz_poly_normalise(R);
}
