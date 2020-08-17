/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void _fmpz_poly_pseudo_div(fmpz * Q, ulong * d, const fmpz * A, slong lenA, 
                          const fmpz * B, slong lenB, const fmpz_preinvn_t inv)
{
    fmpz * R = _fmpz_vec_init(lenA);
    _fmpz_poly_pseudo_divrem(Q, R, d, A, lenA, B, lenB, inv);
    _fmpz_vec_clear(R, lenA);
}

void fmpz_poly_pseudo_div(fmpz_poly_t Q, ulong * d, const fmpz_poly_t A, 
                                                    const fmpz_poly_t B)
{
    slong lenq;
    fmpz * q;

    if (B->length == 0)
    {
        flint_printf("Exception (fmpz_poly_pseudo_div). Division by zero.\n");
        flint_abort();
    }
    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        *d = 0;
        return;
    }

    lenq = A->length - B->length + 1;
    if (Q == A || Q == B)
        q = _fmpz_vec_init(lenq);
    else
    {
        fmpz_poly_fit_length(Q, lenq);
        q = Q->coeffs;
    }

    _fmpz_poly_pseudo_div(q, d, A->coeffs, A->length, 
                                B->coeffs, B->length, NULL);

    if (Q == A || Q == B)
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc  = lenq;
        Q->length = lenq;
    }
    else
        _fmpz_poly_set_length(Q, lenq);
}
