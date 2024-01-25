/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*
   Assumes poly1 and poly2 are not length 0 and 0 < n <= len1 + len2 - 1.
 */
void
_fmpz_poly_mullow_classical(fmpz * res, const fmpz * poly1, slong len1,
                                        const fmpz * poly2, slong len2, slong n)
{
    slong i, top1, top2;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (n == 1)
    {
        fmpz_mul(res, poly1, poly2);
        return;
    }

    if (len1 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly2, n, poly1);
        return;
    }

    if (len2 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly1, n, poly2);
        return;
    }

    fmpz_mul(res, poly1, poly2);

    for (i = 1; i < n; i++)
    {
        top1 = FLINT_MIN(len1 - 1, i);
        top2 = FLINT_MIN(len2 - 1, i);

        _fmpz_vec_dot_general(res + i, NULL, 0, poly1 + i - top2, poly2 + i - top1, 1, top1 + top2 - i + 1);
    }
}

void
fmpz_poly_mullow_classical(fmpz_poly_t res, const fmpz_poly_t poly1,
                                            const fmpz_poly_t poly2, slong n)
{
    slong len_out;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (n > len_out)
        n = len_out;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, n);
        _fmpz_poly_mullow_classical(t->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, n);
        _fmpz_poly_mullow_classical(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, n);
    }

    _fmpz_poly_set_length(res, n);
    _fmpz_poly_normalise(res);
}
