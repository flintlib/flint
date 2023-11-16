/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2012, 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"

int
_gr_poly_compose_series_brent_kung(gr_ptr res, gr_srcptr poly1, slong len1,
                            gr_srcptr poly2, slong len2, slong n, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    gr_mat_t A, B, C;
    gr_ptr t, h;
    slong i, m;

    if (n == 1)
       return gr_set(res, poly1, ctx);

    m = n_sqrt(n) + 1;

    gr_mat_init(A, m, n, ctx);
    gr_mat_init(B, m, m, ctx);
    gr_mat_init(C, m, n, ctx);

    GR_TMP_INIT_VEC(h, 2 * n, ctx);
    t = GR_ENTRY(h, n, sz);

    /* Set rows of B to the segments of poly1 */
    for (i = 0; i < len1 / m; i++)
        status |= _gr_vec_set(B->rows[i], GR_ENTRY(poly1, i * m, sz), m, ctx);
    status |= _gr_vec_set(B->rows[i], GR_ENTRY(poly1, i * m, sz), len1 % m, ctx);

    /* Set rows of A to powers of poly2 */
    status |= gr_one(A->rows[0], ctx);
    status |= _gr_vec_set(A->rows[1], poly2, len2, ctx);

    /* Prefer squaring for powers? todo: the ring should know */
    if (len2 >= n && (gr_ctx_is_finite(ctx) == T_TRUE || gr_ctx_has_real_prec(ctx) == T_TRUE))
    {
        for (i = 2; i < m; i++)
            status |= _gr_poly_mullow(A->rows[i], A->rows[(i + 1) / 2], n, A->rows[i / 2], n, n, ctx);
    }
    else
    {
        for (i = 2; i < m; i++)
            status |= _gr_poly_mullow(A->rows[i], A->rows[i - 1], n, poly2, len2, n, ctx);
    }

    status |= gr_mat_mul(C, B, A, ctx);

    /* Evaluate block composition using the Horner scheme */
    status |= _gr_vec_set(res, C->rows[m - 1], n, ctx);
    status |= _gr_poly_mullow(h, A->rows[m - 1], n, poly2, len2, n, ctx);

    for (i = m - 2; i >= 0; i--)
    {
        status |= _gr_poly_mullow(t, res, n, h, n, n, ctx);
        status |= _gr_poly_add(res, t, n, C->rows[i], n, ctx);
    }

    GR_TMP_CLEAR_VEC(h, 2 * n, ctx);

    gr_mat_clear(A, ctx);
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);

    return status;
}

int
gr_poly_compose_series_brent_kung(gr_poly_t res,
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
        status = _gr_poly_compose_series_brent_kung(res->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
        _gr_poly_set_length(res, lenr, ctx);
        _gr_poly_normalise(res, ctx);
    }
    else
    {
        gr_poly_t t;
        gr_poly_init2(t, lenr, ctx);
        status = _gr_poly_compose_series_brent_kung(t->coeffs, poly1->coeffs, len1,
                                        poly2->coeffs, len2, lenr, ctx);
        _gr_poly_set_length(t, lenr, ctx);
        _gr_poly_normalise(t, ctx);
        gr_poly_swap(res, t, ctx);
        gr_poly_clear(t, ctx);
    }

    return status;
}
