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
_gr_poly_compose_series_horner(gr_ptr res, gr_srcptr poly1, slong len1,
                            gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    slong i = len1 - 1;
    slong lenr;
    slong sz = ctx->sizeof_elem;
    gr_ptr t;
    int status = GR_SUCCESS;

    if (n == 1)
        return gr_set(res, poly1, ctx);

    lenr = len2;
    status |= _gr_vec_mul_scalar(res, poly2, len2, GR_ENTRY(poly1, i, sz), ctx);
    i--;
    status |= gr_add(res, res, GR_ENTRY(poly1, i, sz), ctx);

    if (i > 0)
    {
        GR_TMP_INIT_VEC(t, n, ctx);

        while (i > 0)
        {
            i--;
            if (lenr + len2 - 1 < n)
            {
                status |= _gr_poly_mul(t, res, lenr, poly2, len2, ctx);
                lenr = lenr + len2 - 1;
            }
            else
            {
                status |= _gr_poly_mullow(t, res, lenr, poly2, len2, n, ctx);
                lenr = n;
            }
            status |= _gr_poly_add(res, t, lenr, GR_ENTRY(poly1, i, sz), 1, ctx);
        }

        GR_TMP_CLEAR_VEC(t, n, ctx);
    }

    status |= _gr_vec_zero(GR_ENTRY(res, lenr, sz), n - lenr, ctx);

    return status;
}

int
gr_poly_compose_series_horner(gr_poly_t res,
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
        status = _gr_poly_compose_series_horner(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
        _gr_poly_set_length(res, lenr, ctx);
        _gr_poly_normalise(res, ctx);
    }
    else
    {
        gr_poly_t t;
        gr_poly_init2(t, lenr, ctx);
        status = _gr_poly_compose_series_horner(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
        _gr_poly_set_length(t, lenr, ctx);
        _gr_poly_normalise(t, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }

    return status;
}
