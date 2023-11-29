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

void _fmpz_poly_pseudo_rem(fmpz * R, ulong * d, const fmpz * A, slong lenA,
                          const fmpz * B, slong lenB, const fmpz_preinvn_t inv)
{
    fmpz * Q = _fmpz_vec_init(lenA + lenB - 1);
    _fmpz_poly_pseudo_divrem(Q, R, d, A, lenA, B, lenB, inv);
    _fmpz_vec_clear(Q, lenA + lenB - 1);
}

void fmpz_poly_pseudo_rem(fmpz_poly_t R, ulong * d, const fmpz_poly_t A,
                                                    const fmpz_poly_t B)
{
    slong lenr;
    fmpz * r;

    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_pseudo_rem). Division by zero.\n");
    }
    if (A->length < B->length)
    {
        fmpz_poly_set(R, A);
        *d = 0;
        return;
    }

    lenr = A->length;
    if (R == B || R == A)
        r = _fmpz_vec_init(lenr);
    else
    {
        fmpz_poly_fit_length(R, lenr);
        r = R->coeffs;
    }

    _fmpz_poly_pseudo_rem(r, d, A->coeffs, A->length,
                                B->coeffs, B->length, NULL);

    for (lenr = B->length - 2; (lenr >= 0) && !r[lenr]; lenr--) ;
    lenr++;

    if (R == B || R == A)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc  = A->length;
        R->length = lenr;
    }
    else
        _fmpz_poly_set_length(R, lenr);
}
