/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"
#include "acb_mat.h"

void
_acb_poly_compose_series_brent_kung(acb_ptr res,
    acb_srcptr poly1, slong len1,
    acb_srcptr poly2, slong len2, slong n, slong prec)
{
    acb_mat_t A, B, C;
    acb_ptr t, h;
    slong i, m;

    if (n == 1)
    {
        acb_set(res, poly1);
        return;
    }

    m = n_sqrt(n) + 1;

    acb_mat_init(A, m, n);
    acb_mat_init(B, m, m);
    acb_mat_init(C, m, n);

    h = _acb_vec_init(n);
    t = _acb_vec_init(n);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        _acb_vec_set(B->rows[i], poly1 + i*m, m);
    _acb_vec_set(B->rows[i], poly1 + i*m, len1 % m);

    /* Set rows of A to powers of poly2 */
    acb_set_ui(A->rows[0] + 0, UWORD(1));
    _acb_vec_set(A->rows[1], poly2, len2);
    for (i = 2; i < m; i++)
        _acb_poly_mullow(A->rows[i], A->rows[(i + 1) / 2], n, A->rows[i / 2], n, n, prec);

    acb_mat_mul(C, B, A, prec);

    /* Evaluate block composition using the Horner scheme */
    _acb_vec_set(res, C->rows[m - 1], n);
    _acb_poly_mullow(h, A->rows[m - 1], n, poly2, len2, n, prec);

    for (i = m - 2; i >= 0; i--)
    {
        _acb_poly_mullow(t, res, n, h, n, n, prec);
        _acb_poly_add(res, t, n, C->rows[i], n, prec);
    }

    _acb_vec_clear(h, n);
    _acb_vec_clear(t, n);

    acb_mat_clear(A);
    acb_mat_clear(B);
    acb_mat_clear(C);
}

void
acb_poly_compose_series_brent_kung(acb_poly_t res,
                    const acb_poly_t poly1,
                    const acb_poly_t poly2, slong n, slong prec)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;

    if (len2 != 0 && !acb_is_zero(poly2->coeffs))
    {
        flint_printf("exception: compose_series: inner "
                "polynomial must have zero constant term\n");
        flint_abort();
    }

    if (len1 == 0 || n == 0)
    {
        acb_poly_zero(res);
        return;
    }

    if (len2 == 0 || len1 == 1)
    {
        acb_poly_set_acb(res, poly1->coeffs);
        return;
    }

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        acb_poly_fit_length(res, lenr);
        _acb_poly_compose_series_brent_kung(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _acb_poly_set_length(res, lenr);
        _acb_poly_normalise(res);
    }
    else
    {
        acb_poly_t t;
        acb_poly_init2(t, lenr);
        _acb_poly_compose_series_brent_kung(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _acb_poly_set_length(t, lenr);
        _acb_poly_normalise(t);
        acb_poly_swap(res, t);
        acb_poly_clear(t);
    }
}
