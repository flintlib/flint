/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, slong len1,
                            arb_srcptr poly2, slong len2, slong n, slong prec)
{
    if (len2 == 1)
    {
        arb_set_round(res, poly1, prec);
        _arb_vec_zero(res + 1, n - 1);
    }
    else if (_arb_vec_is_zero(poly2 + 1, len2 - 2))  /* poly2 is a monomial */
    {
        slong i, j;
        arb_t t;

        arb_init(t);
        arb_set(t, poly2 + len2 - 1);
        arb_set_round(res, poly1, prec);

        for (i = 1, j = len2 - 1; i < len1 && j < n; i++, j += len2 - 1)
        {
            arb_mul(res + j, poly1 + i, t, prec);

            if (i + 1 < len1 && j + len2 - 1 < n)
                arb_mul(t, t, poly2 + len2 - 1, prec);
        }

        if (len2 != 2)
            for (i = 1; i < n; i++)
                if (i % (len2 - 1) != 0)
                    arb_zero(res + i);

        arb_clear(t);
    }
    else if (len1 < 6 || n < 6)
    {
        _arb_poly_compose_series_horner(res, poly1, len1, poly2, len2, n, prec);
    }
    else
    {
        _arb_poly_compose_series_brent_kung(res, poly1, len1, poly2, len2, n, prec);
    }
}

void
arb_poly_compose_series(arb_poly_t res,
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
        _arb_poly_compose_series(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _arb_poly_set_length(res, lenr);
        _arb_poly_normalise(res);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init2(t, lenr);
        _arb_poly_compose_series(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _arb_poly_set_length(t, lenr);
        _arb_poly_normalise(t);
        arb_poly_swap(res, t);
        arb_poly_clear(t);
    }
}
