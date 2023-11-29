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
_fmpz_poly_div_basecase(fmpz * Q, fmpz * R, const fmpz * A, slong lenA,
                        const fmpz * B, slong lenB, int exact)
{
    const fmpz * leadB = B + lenB - 1;
    slong B1, iQ = lenA - lenB;
    slong alloc;
    fmpz_t r;
    int res = 1;

    while (lenA >= lenB && fmpz_cmpabs(A + lenA - 1, leadB) < 0)
    {
        if (exact && !fmpz_is_zero(A + lenA - 1))
           return 0;

        fmpz_zero(Q + iQ);
        iQ--;
        lenA--;
    }

    if (lenA < lenB)
        return 1;

    alloc = (R == NULL) ? lenA : 0;
    if (alloc)
        R = _fmpz_vec_init(alloc);
    if (R != A)
        _fmpz_vec_set(R + lenB - 1, A + lenB - 1, lenA - lenB + 1);

    B1 = lenB - 1;

    if (exact)
        fmpz_init(r);

    while (lenA >= lenB)
    {
        if (fmpz_cmpabs(R + lenA - 1, leadB) < 0)
        {
            if (exact && !fmpz_is_zero(R + lenA - 1))
            {
                res = 0;
                goto cleanup;
            }

            fmpz_zero(Q + iQ);
        } else
        {
            if (exact)
            {
                fmpz_fdiv_qr(Q + iQ, r, R + lenA - 1, leadB);

                if (!fmpz_is_zero(r))
                {
                    res = 0;
                    goto cleanup;
                }
            } else
                fmpz_fdiv_q(Q + iQ, R + lenA - 1, leadB);

            _fmpz_vec_scalar_submul_fmpz(R + lenA - B1 - 1, B, B1, Q + iQ);
        }

        if (B1 >= lenA - lenB + 1)
        {
            B++;
            B1--;
        }

        lenA--;
        iQ--;
    }

cleanup:

    if (exact)
        fmpz_clear(r);

    if (alloc)
        _fmpz_vec_clear(R, alloc);

    return res;
}

void
fmpz_poly_div_basecase(fmpz_poly_t Q,
                       const fmpz_poly_t A, const fmpz_poly_t B)
{
    slong lenq;
    fmpz *q;

    if (B->length == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_div_basecase). Division by zero.\n");
    }
    if (A->length < B->length)
    {
        fmpz_poly_zero(Q);
        return;
    }

    lenq = A->length - B->length + 1;
    if ((Q == A) || (Q == B))
        q = _fmpz_vec_init(lenq);
    else
    {
        fmpz_poly_fit_length(Q, lenq);
        q = Q->coeffs;
    }

    _fmpz_poly_div_basecase(q, NULL, A->coeffs, A->length,
                                     B->coeffs, B->length, 0);

    if ((Q == A) || (Q == B))
    {
        _fmpz_vec_clear(Q->coeffs, Q->alloc);
        Q->coeffs = q;
        Q->alloc = lenq;
        Q->length = lenq;
    }
    else
        _fmpz_poly_set_length(Q, lenq);

    _fmpz_poly_normalise(Q);
}
