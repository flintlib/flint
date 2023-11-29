/*
    Copyright (C) 2012, 2013 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "gr_poly.h"

void
_arb_poly_compose_series(arb_ptr res, arb_srcptr poly1, slong len1,
                            arb_srcptr poly2, slong len2, slong n, slong prec)
{
    if (len2 == 1)
    {
        arb_set_round(res, poly1, prec);
        _arb_vec_zero(res + 1, n - 1);
    }
    else if (!_arb_vec_is_finite(poly1, len1) || !_arb_vec_is_finite(poly2, len2))
    {
        /* find k such that the first k coefficients of both poly1 and
         * poly2 are finite */
        slong k = 0;

        while (((k >= len1) || arb_is_finite(poly1 + k)) && ((k >= len2) || arb_is_finite(poly2 + k)))
            k += 1;

        if (k > 0)
        {
            gr_ctx_t ctx;
            gr_ctx_init_real_arb(ctx, prec);
            GR_MUST_SUCCEED(_gr_poly_compose_series(res, poly1, FLINT_MIN(len1, k), poly2, FLINT_MIN(len2, k), FLINT_MIN(n, k), ctx));
            _arb_vec_indeterminate(res + k, n - k);
        }
        else
        {
            _arb_vec_indeterminate(res, n);
        }
    }
    else
    {
        gr_ctx_t ctx;
        gr_ctx_init_real_arb(ctx, prec);
        GR_MUST_SUCCEED(_gr_poly_compose_series(res, poly1, len1, poly2, len2, n, ctx));
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
        flint_throw(FLINT_ERROR, "arb_poly_compose_series: "
                "inner polynomial must have zero constant term\n");
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
