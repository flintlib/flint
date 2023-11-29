/*
    Copyright (C) 2008, 2009 William Hart
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

int
_fmpz_poly_divrem_basecase(fmpz * Q, fmpz * R, const fmpz * A, slong lenA,
                           const fmpz * B, slong lenB, int exact)
{
    const fmpz * leadB = B + (lenB - 1);
    slong iQ, iR;
    fmpz_t r;

    if (exact)
       fmpz_init(r);

    if (R != A)
        _fmpz_vec_set(R, A, lenA);

    for (iQ = lenA - lenB, iR = lenA - 1; iQ >= 0; iQ--, iR--)
    {
        if (fmpz_cmpabs(R + iR, leadB) < 0)
        {
            if (exact && !fmpz_is_zero(R + iR))
            {
                fmpz_clear(r);
                return 0;
            }

            fmpz_zero(Q + iQ);
        } else
        {
            if (exact)
            {
                fmpz_fdiv_qr(Q + iQ, r, R + iR, leadB);

                if (!fmpz_is_zero(r))
                {
                    fmpz_clear(r);
                    return 0;
                }
            } else
                fmpz_fdiv_q(Q + iQ, R + iR, leadB);

            _fmpz_vec_scalar_submul_fmpz(R + iQ, B, lenB, Q + iQ);
        }
    }

    if (exact)
       fmpz_clear(r);

    return 1;
}

void
fmpz_poly_divrem_basecase(fmpz_poly_t Q, fmpz_poly_t R,
                          const fmpz_poly_t A, const fmpz_poly_t B)
{
    slong lenq, lenr;
    fmpz *q, *r;

    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_divrem_basecase). Division by zero.\n");
    }
    if (Q == R)
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_divrem_basecase): Output arguments Q and R may not be aliased.\n");
    }
    if (A->length < B->length)
    {
        fmpz_poly_set(R, A);
        fmpz_poly_zero(Q);
        return;
    }

    lenq = A->length - B->length + 1;
    lenr = A->length;
    if (Q == A || Q == B)
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

    _fmpz_poly_divrem_basecase(q, r, A->coeffs, A->length,
                                     B->coeffs, B->length, 0);

    if (Q == A || Q == B)
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
        R->alloc = lenr;
        R->length = lenr;
    }
    else
        _fmpz_poly_set_length(R, lenr);

    _fmpz_poly_normalise(Q);
    _fmpz_poly_normalise(R);
}
