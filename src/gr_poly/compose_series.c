/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_compose_series(gr_ptr res, gr_srcptr poly1, slong len1,
                            gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (n == 1)
    {
        status |= gr_set(res, poly1, ctx);
    }
    else if (n == 2)
    {
        status |= gr_set(res, poly1, ctx);
        status |= gr_mul(GR_ENTRY(res, 1, sz), GR_ENTRY(poly1, 1, sz), GR_ENTRY(poly2, 1, sz), ctx);
    }
    else if (len1 == 1 || len2 == 1)
    {
        status |= gr_set(res, poly1, ctx);
        status |= _gr_vec_zero(GR_ENTRY(res, 1, sz), n - 1, ctx);
    }
    else if (n == 3)
    {
        status |= gr_set(res, poly1, ctx);
        status |= gr_mul(GR_ENTRY(res, 1, sz), GR_ENTRY(poly1, 1, sz), GR_ENTRY(poly2, 1, sz), ctx);

        if (len1 == 2 && len2 == 3)
        {
            status |= gr_mul(GR_ENTRY(res, 2, sz), GR_ENTRY(poly1, 1, sz), GR_ENTRY(poly2, 2, sz), ctx);
        }
        else
        {
            status |= gr_sqr(GR_ENTRY(res, 2, sz), GR_ENTRY(poly2, 1, sz), ctx);
            status |= gr_mul(GR_ENTRY(res, 2, sz), GR_ENTRY(res, 2, sz), GR_ENTRY(poly1, 2, sz), ctx);
            if (len2 == 3)
                status |= gr_addmul(GR_ENTRY(res, 2, sz), GR_ENTRY(poly1, 1, sz), GR_ENTRY(poly2, 2, sz), ctx);
        }
    }
    else if (_gr_vec_is_zero(GR_ENTRY(poly2, 1, sz), len2 - 2, ctx) == T_TRUE)  /* poly2 is a monomial */
    {
        slong i, j;
        gr_ptr t;

        GR_TMP_INIT(t, ctx);
        status |= gr_set(t, GR_ENTRY(poly2, len2 - 1, sz), ctx);
        status |= gr_set(res, poly1, ctx);

        for (i = 1, j = len2 - 1; i < len1 && j < n; i++, j += len2 - 1)
        {
            status |= gr_mul(GR_ENTRY(res, j, sz), GR_ENTRY(poly1, i, sz), t, ctx);

            if (i + 1 < len1 && j + len2 - 1 < n)
                status |= gr_mul(t, t, GR_ENTRY(poly2, len2 - 1, sz), ctx);
        }

        if (len2 != 2)
            for (i = 1; i < n; i++)
                if (i % (len2 - 1) != 0)
                    status |= gr_zero(GR_ENTRY(res, i, sz), ctx);

        GR_TMP_CLEAR(t, ctx);
    }
    else if (len1 <= 7 || n <= 7)
    {
        status |= _gr_poly_compose_series_horner(res, poly1, len1, poly2, len2, n, ctx);
    }
    else if (len1 * len1 < n || (len1 - 1) * (len2 - 1) + 1 < 4 * n)  /* todo: tune these cutoffs */
    {
        status |= _gr_poly_compose_series_divconquer(res, poly1, len1, poly2, len2, n, ctx);
    }
    else
    {
        status |= _gr_poly_compose_series_brent_kung(res, poly1, len1, poly2, len2, n, ctx);
    }

    return status;
}

int
gr_poly_compose_series(gr_poly_t res,
                    const gr_poly_t poly1,
                    const gr_poly_t poly2, slong n, gr_ctx_t ctx)
{
    slong len1 = poly1->length;
    slong len2 = poly2->length;
    slong lenr;
    int status;

    if (len2 != 0)
    {
        truth_t is_zero = gr_is_zero(poly2->coeffs, ctx);

        if (is_zero == T_FALSE)
            return GR_DOMAIN;
        if (is_zero == T_UNKNOWN)
            return GR_UNABLE;
    }

    if (len1 == 0 || n == 0)
        return gr_poly_zero(res, ctx);

    if (len2 == 0 || len1 == 1)
        return gr_poly_set_scalar(res, poly1->coeffs, ctx);

    lenr = FLINT_MIN((len1 - 1) * (len2 - 1) + 1, n);
    len1 = FLINT_MIN(len1, lenr);
    len2 = FLINT_MIN(len2, lenr);

    if ((res != poly1) && (res != poly2))
    {
        gr_poly_fit_length(res, lenr, ctx);
        status = _gr_poly_compose_series(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
        _gr_poly_set_length(res, lenr, ctx);
        _gr_poly_normalise(res, ctx);
    }
    else
    {
        gr_poly_t t;
        gr_poly_init2(t, lenr, ctx);
        status = _gr_poly_compose_series(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
        _gr_poly_set_length(t, lenr, ctx);
        _gr_poly_normalise(t, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }

    return status;
}
