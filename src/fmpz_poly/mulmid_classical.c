/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
_fmpz_poly_mulmid_classical(fmpz * res, const fmpz * poly1, slong len1,
                                        const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    slong i, top1, top2, start, stop, len;

    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nhi != 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

    len = nhi - nlo;

    if (len1 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly2 + nlo, len, poly1);
        return;
    }

    if (len2 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly1 + nlo, len, poly2);
        return;
    }

    if (poly1 == poly2 && len1 == len2)
    {
        if (nlo == 0)
            fmpz_mul(res, poly1, poly1);

        for (i = FLINT_MAX(nlo, 1); i < FLINT_MIN(nhi, 2 * len1 - 2); i++)
        {
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            _fmpz_vec_dot_general(res + i - nlo, NULL, 0,
                poly1 + start, poly1 + i - stop, 1, stop - start + 1);
            fmpz_mul_2exp(res + i - nlo, res + i - nlo, 1);

            if (i % 2 == 0)
                fmpz_addmul(res + i - nlo, poly1 + i / 2, poly1 + i / 2);
        }

        if (nhi >= 2 * len1 - 1)
            fmpz_mul(res + 2 * len1 - 2 - nlo, poly1 + len1 - 1, poly1 + len1 - 1);
    }
    else
    {
        for (i = nlo; i < nhi; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);

            _fmpz_vec_dot_general(res + i - nlo, NULL, 0,
                poly1 + i - top2, poly2 + i - top1, 1, top1 + top2 - i + 1);
        }
    }
}

void
fmpz_poly_mulmid_classical(fmpz_poly_t res, const fmpz_poly_t poly1,
                                            const fmpz_poly_t poly2, slong nlo, slong nhi)
{
    slong len1 = fmpz_poly_length(poly1);
    slong len2 = fmpz_poly_length(poly2);
    slong len;

    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi >= 0);

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
    {
        fmpz_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, len);
        _fmpz_poly_mulmid_classical(t->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, nlo, nhi);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, len);
        _fmpz_poly_mulmid_classical(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, nlo, nhi);
    }

    _fmpz_poly_set_length(res, len);
    _fmpz_poly_normalise(res);
}
