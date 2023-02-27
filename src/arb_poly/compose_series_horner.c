/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"

void
_arb_poly_compose_series_horner(arb_ptr res, arb_srcptr poly1, slong len1,
                            arb_srcptr poly2, slong len2, slong n, slong prec)
{
    if (n == 1)
    {
        arb_set(res, poly1);
    }
    else
    {
        slong i = len1 - 1;
        slong lenr;

        arb_ptr t = _arb_vec_init(n);

        lenr = len2;
        _arb_vec_scalar_mul(res, poly2, len2, poly1 + i, prec);
        i--;
        arb_add(res, res, poly1 + i, prec);

        while (i > 0)
        {
            i--;
            if (lenr + len2 - 1 < n)
            {
                _arb_poly_mul(t, res, lenr, poly2, len2, prec);
                lenr = lenr + len2 - 1;
            }
            else
            {
                _arb_poly_mullow(t, res, lenr, poly2, len2, n, prec);
                lenr = n;
            }
            _arb_poly_add(res, t, lenr, poly1 + i, 1, prec);
        }

        _arb_vec_zero(res + lenr, n - lenr);
        _arb_vec_clear(t, n);
    }
}

void
arb_poly_compose_series_horner(arb_poly_t res,
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
        _arb_poly_compose_series_horner(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _arb_poly_set_length(res, lenr);
        _arb_poly_normalise(res);
    }
    else
    {
        arb_poly_t t;
        arb_poly_init2(t, lenr);
        _arb_poly_compose_series_horner(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, prec);
        _arb_poly_set_length(t, lenr);
        _arb_poly_normalise(t);
        arb_poly_swap(res, t);
        arb_poly_clear(t);
    }
}
