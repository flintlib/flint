/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "acb_poly.h"

void
_acb_poly_compose_series(acb_ptr res, acb_srcptr poly1, slong len1,
                            acb_srcptr poly2, slong len2, slong n, slong prec)
{
    if (len2 == 1)
    {
        acb_set_round(res, poly1, prec);
        _acb_vec_zero(res + 1, n - 1);
    }
    else if (_acb_vec_is_zero(poly2 + 1, len2 - 2))  /* poly2 is a monomial */
    {
        slong i, j;
        acb_t t;

        acb_init(t);
        acb_set(t, poly2 + len2 - 1);
        acb_set_round(res, poly1, prec);

        for (i = 1, j = len2 - 1; i < len1 && j < n; i++, j += len2 - 1)
        {
            acb_mul(res + j, poly1 + i, t, prec);

            if (i + 1 < len1 && j + len2 - 1 < n)
                acb_mul(t, t, poly2 + len2 - 1, prec);
        }

        if (len2 != 2)
            for (i = 1; i < n; i++)
                if (i % (len2 - 1) != 0)
                    acb_zero(res + i);

        acb_clear(t);
    }
    else if (len1 < 6 || n < 6)
    {
        _acb_poly_compose_series_horner(res, poly1, len1, poly2, len2, n, prec);
    }
    else
    {
        _acb_poly_compose_series_brent_kung(res, poly1, len1, poly2, len2, n, prec);
    }
}

void
acb_poly_compose_series(acb_poly_t res,
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
        _acb_poly_compose_series(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _acb_poly_set_length(res, lenr);
        _acb_poly_normalise(res);
    }
    else
    {
        acb_poly_t t;
        acb_poly_init2(t, lenr);
        _acb_poly_compose_series(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _acb_poly_set_length(t, lenr);
        _acb_poly_normalise(t);
        acb_poly_swap(res, t);
        acb_poly_clear(t);
    }
}
