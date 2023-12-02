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

void
_fmpz_poly_pseudo_divrem_cohen(fmpz * Q, fmpz * R, const fmpz * A,
                               slong lenA, const fmpz * B, slong lenB)
{
    const fmpz * leadB = B + (lenB - 1);
    slong e, lenQ;
    fmpz_t pow;

    if (lenB == 1)
    {
        fmpz_init(pow);
        fmpz_pow_ui(pow, leadB, lenA - 1);
        _fmpz_vec_scalar_mul_fmpz(Q, A, lenA, pow);
        _fmpz_vec_zero(R, lenA);
        fmpz_clear(pow);
        return;
    }

    lenQ = lenA - lenB + 1;
    _fmpz_vec_zero(Q, lenQ);
    if (R != A)
        _fmpz_vec_set(R, A, lenA);
    e = lenA - lenB;

    /* Unroll the first run of the while loop */
    {
        fmpz_set(Q + (lenQ - 1), R + (lenA - 1));

        _fmpz_vec_scalar_mul_fmpz(R, R, lenA - 1, leadB);
        _fmpz_vec_scalar_submul_fmpz(R + (lenA - lenB), B, lenB - 1, R + (lenA - 1));
        fmpz_zero(R + (lenA - 1));

        lenA--;
        FMPZ_VEC_NORM(R, lenA);
    }
    while (lenA >= lenB)
    {
        _fmpz_vec_scalar_mul_fmpz(Q, Q, lenQ, leadB);
        fmpz_add(Q + (lenA - lenB), Q + (lenA - lenB), R + (lenA - 1));

        _fmpz_vec_scalar_mul_fmpz(R, R, lenA - 1, leadB);
        _fmpz_vec_scalar_submul_fmpz(R + lenA - lenB, B, lenB - 1, R + (lenA - 1));
        fmpz_zero(R + (lenA - 1));

        lenA--;
        FMPZ_VEC_NORM(R, lenA);

        e--;
    }

    fmpz_init(pow);
    fmpz_pow_ui(pow, leadB, e);
    _fmpz_vec_scalar_mul_fmpz(Q, Q, lenQ, pow);
    _fmpz_vec_scalar_mul_fmpz(R, R, lenA, pow);
    fmpz_clear(pow);
}

void
fmpz_poly_pseudo_divrem_cohen(fmpz_poly_t Q, fmpz_poly_t R,
                              const fmpz_poly_t A, const fmpz_poly_t B)
{
    slong lenq, lenr;
    fmpz *q, *r;

    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_pseudo_divrem_cohen): Division by zero.\n");
    }
    if (Q == R)
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_pseudo_divrem_cohen): "
                "Output arguments Q and R may not be aliased.\n");
    }
    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        fmpz_poly_set(R, A);
        return;
    }

    lenq = A->length - B->length + 1;
    lenr = A->length;
    if ((Q == A) || (Q == B))
        q = _fmpz_vec_init(lenq);
    else
    {
        fmpz_poly_fit_length(Q, lenq);
        q = Q->coeffs;
    }
    if (R == B)
        r = _fmpz_vec_init(lenr);
    else
    {
        fmpz_poly_fit_length(R, lenr);
        r = R->coeffs;
    }

    _fmpz_poly_pseudo_divrem_cohen(q, r, A->coeffs, A->length, B->coeffs, B->length);

    for (lenr = B->length - 1; (lenr >= 0) && r[lenr] == WORD(0); lenr--) ;
    lenr++;

    if ((Q == A) || (Q == B))
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc = lenq;
        Q->length = lenq;
    }
    else
        _fmpz_poly_set_length(Q, lenq);
    if (R == B)
    {
        _fmpz_vec_clear(R->coeffs, R->alloc);
        R->coeffs = r;
        R->alloc = A->length;
        R->length = lenr;
    }
    else
        _fmpz_poly_set_length(R, lenr);
}

