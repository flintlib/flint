/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_poly.h"

int
_gr_poly_mullow_generic(gr_ptr res,
    gr_srcptr poly1, slong len1,
    gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    int status;
    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (n == 1)
        return gr_mul(res, poly1, poly2, ctx);

    if (len1 == 1)
        return _gr_vec_mul_scalar(res, poly2, n, poly1, ctx);

    if (len2 == 1)
        return _gr_vec_mul_scalar(res, poly1, n, poly2, ctx);

    /* Squaring */
    if (poly1 == poly2 && len1 == len2)
    {
        slong i, start, stop, sz;
        sz = ctx->sizeof_elem;

        status = GR_SUCCESS;

        /* todo: double, square, addmul */

        status |= gr_sqr(res, poly1, ctx);
        status |= gr_mul(GR_ENTRY(res, 1, sz), poly1, GR_ENTRY(poly1, 1, sz), ctx);
        status |= gr_mul_two(GR_ENTRY(res, 1, sz), GR_ENTRY(res, 1, sz), ctx);

        for (i = 2; i < FLINT_MIN(n, 2 * len1 - 3); i++)
        {
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            status |= _gr_vec_dot_rev(GR_ENTRY(res, i, sz), NULL, 0, GR_ENTRY(poly1, start, sz), GR_ENTRY(poly1, i - stop, sz), stop - start + 1, ctx);
            status |= gr_mul_two(GR_ENTRY(res, i, sz), GR_ENTRY(res, i, sz), ctx);

            if (i % 2 == 0 && i / 2 < len1)   /* should be addsqr? */
                status |= gr_addmul(GR_ENTRY(res, i, sz), GR_ENTRY(poly1, i / 2, sz), GR_ENTRY(poly1, i / 2, sz), ctx);
        }

        if (len1 > 2 && n >= 2 * len1 - 2)
        {
            status |= gr_mul(GR_ENTRY(res, 2 * len1 - 3, sz), GR_ENTRY(poly1, len1 - 1, sz), GR_ENTRY(poly1, len1 - 2, sz), ctx);
            status |= gr_mul_two(GR_ENTRY(res, 2 * len1 - 3, sz), GR_ENTRY(res, 2 * len1 - 3, sz), ctx);
        }

        if (n >= 2 * len1 - 1)
            status |= gr_sqr(GR_ENTRY(res, 2 * len1 - 2, sz), GR_ENTRY(poly1, len1 - 1, sz), ctx);

        return status;
    }
    else
    {
        slong i, top1, top2, sz;
        sz = ctx->sizeof_elem;

        status = gr_mul(res, poly1, poly2, ctx);

        for (i = 1; i < n; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);

            status |= _gr_vec_dot_rev(GR_ENTRY(res, i, sz), NULL, 0, GR_ENTRY(poly1, i - top2, sz), GR_ENTRY(poly2, i - top1, sz), top1 + top2 - i + 1, ctx);
        }

        return status;
    }
}

int
gr_poly_mullow(gr_poly_t res, const gr_poly_t poly1,
                                            const gr_poly_t poly2,
                                                slong n, gr_ctx_t ctx)
{
    slong len_out;
    int status;

    if (poly1->length == 0 || poly2->length == 0 || n == 0)
        return gr_poly_zero(res, ctx);

    len_out = poly1->length + poly2->length - 1;
    n = FLINT_MIN(n, len_out);

    if (res == poly1 || res == poly2)
    {
        gr_poly_t t;
        gr_poly_init2(t, n, ctx);
        status = _gr_poly_mullow(t->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }
    else
    {
        gr_poly_fit_length(res, n, ctx);
        status = _gr_poly_mullow(res->coeffs, poly1->coeffs, poly1->length, poly2->coeffs, poly2->length, n, ctx);
    }

    _gr_poly_set_length(res, n, ctx);
    _gr_poly_normalise(res, ctx);
    return status;
}
