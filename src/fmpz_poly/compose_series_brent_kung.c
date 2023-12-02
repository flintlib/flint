/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"

void
_fmpz_poly_compose_series_brent_kung(fmpz * res, const fmpz * poly1, slong len1,
                                      const fmpz * poly2, slong len2, slong n)
{
    fmpz_mat_t A, B, C;
    fmpz *t, *h;
    slong i, m;

    if (n == 1)
    {
        fmpz_set(res, poly1);
        return;
    }

    m = n_sqrt(n) + 1;

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(B, m, m);
    fmpz_mat_init(C, m, n);

    h = _fmpz_vec_init(n);
    t = _fmpz_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _fmpz_vec_set(B->rows[i], poly1 + i*m, m);
    _fmpz_vec_set(B->rows[i], poly1 + i*m, len1 % m);

    /* Set rows of A to powers of poly2 */
    fmpz_one(A->rows[0]);
    _fmpz_vec_set(A->rows[1], poly2, len2);
    for (i = 2; i < m; i++)
        _fmpz_poly_mullow(A->rows[i], A->rows[i-1], n, poly2, len2, n);

    fmpz_mat_mul(C, B, A);

    /* Evaluate block composition using the Horner scheme */
    _fmpz_vec_set(res, C->rows[m - 1], n);
    _fmpz_poly_mullow(h, A->rows[m - 1], n, poly2, len2, n);

    for (i = m - 2; i >= 0; i--)
    {
        _fmpz_poly_mullow(t, res, n, h, n, n);
        _fmpz_poly_add(res, t, n, C->rows[i], n);
    }

    _fmpz_vec_clear(h, n);
    _fmpz_vec_clear(t, n);

    fmpz_mat_clear(A);
    fmpz_mat_clear(B);
    fmpz_mat_clear(C);
}

void
fmpz_poly_compose_series_brent_kung(fmpz_poly_t res,
                    const fmpz_poly_t poly1, const fmpz_poly_t poly2, slong n)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;

    if (len2 != 0 && !fmpz_is_zero(poly2->coeffs))
    {
        flint_throw(FLINT_ERROR, "(fmpz_poly_compose_series_brent_kung): "
                "Inner polynomial must have zero constant term.\n");
    }

    if (len1 == 0 || n == 0)
    {
        fmpz_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        fmpz_poly_set_fmpz(res, poly1->coeffs);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        fmpz_poly_fit_length(res, lenr);
        _fmpz_poly_compose_series_brent_kung(res->coeffs, poly1->coeffs, len1,
                                               poly2->coeffs, len2, lenr);
        _fmpz_poly_set_length(res, lenr);
        _fmpz_poly_normalise(res);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, lenr);
        _fmpz_poly_compose_series_brent_kung(t->coeffs, poly1->coeffs, len1,
                                             poly2->coeffs, len2, lenr);
        _fmpz_poly_set_length(t, lenr);
        _fmpz_poly_normalise(t);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
