/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "arb_mat.h"

void
_arb_poly_compose_series_brent_kung(arb_ptr res,
    arb_srcptr poly1, slong len1,
    arb_srcptr poly2, slong len2, slong n, slong prec)
{
    arb_mat_t A, B, C;
    arb_ptr t, h;
    slong i, m;

    if (n == 1)
    {
        arb_set(res, poly1);
        return;
    }

    m = n_sqrt(n) + 1;

    arb_mat_init(A, m, n);
    arb_mat_init(B, m, m);
    arb_mat_init(C, m, n);

    h = _arb_vec_init(n);
    t = _arb_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _arb_vec_set(B->rows[i], poly1 + i*m, m);
    _arb_vec_set(B->rows[i], poly1 + i*m, len1 % m);

    /* Set rows of A to powers of poly2 */
    arb_set_ui(A->rows[0] + 0, UWORD(1));
    _arb_vec_set(A->rows[1], poly2, len2);
    for (i = 2; i < m; i++)
        _arb_poly_mullow(A->rows[i], A->rows[(i + 1) / 2], n, A->rows[i / 2], n, n, prec);

    arb_mat_mul(C, B, A, prec);

    /* Evaluate block composition using the Horner scheme */
    _arb_vec_set(res, C->rows[m - 1], n);
    _arb_poly_mullow(h, A->rows[m - 1], n, poly2, len2, n, prec);

    for (i = m - 2; i >= 0; i--)
    {
        _arb_poly_mullow(t, res, n, h, n, n, prec);
        _arb_poly_add(res, t, n, C->rows[i], n, prec);
    }

    _arb_vec_clear(h, n);
    _arb_vec_clear(t, n);

    arb_mat_clear(A);
    arb_mat_clear(B);
    arb_mat_clear(C);
}

void
arb_poly_compose_series_brent_kung(arb_poly_t res,
                    const arb_poly_t poly1,
                    const arb_poly_t poly2, slong n, slong prec)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;

    if (len2 != 0 && !arb_is_zero(poly2->coeffs))
    {
        flint_printf("exception: compose_series: inner "
                "polynomial must have zero constant term\n");
        flint_abort();
    }

    if (len1 == 0 || n == 0)
    {
        arb_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        arb_poly_set_arb(res, poly1->coeffs);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        arb_poly_fit_length(res, lenr);
        _arb_poly_compose_series_brent_kung(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _arb_poly_set_length(res, lenr);
        _arb_poly_normalise(res);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init2(t, lenr);
        _arb_poly_compose_series_brent_kung(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _arb_poly_set_length(t, lenr);
        _arb_poly_normalise(t);
        arb_poly_swap(res, t);
        arb_poly_clear(t);
    }
}
