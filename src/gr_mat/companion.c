/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_poly.h"

int
_gr_mat_companion(gr_mat_t res, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, n;
    gr_ptr c;
    slong sz = ctx->sizeof_elem;

    n = gr_mat_nrows(res, ctx);

    if (n == 0)
        return status;

    for (i = 0; i < n - 1; i++)
        for (j = 0; j < n; j++)
            status |= gr_set_ui(gr_mat_entry_ptr(res, i, j, ctx), (i + 1) == j, ctx);

    GR_TMP_INIT(c, ctx);
    status |= gr_inv(c, GR_ENTRY(poly, n, sz), ctx);
    status |= gr_neg(c, c, ctx);
    for (j = 0; j < n; j++)
        status |= gr_mul(gr_mat_entry_ptr(res, n - 1, j, ctx), GR_ENTRY(poly, j, sz), c, ctx);
    GR_TMP_CLEAR(c, ctx);

    return status;
}


int
gr_mat_companion(gr_mat_t res, const gr_poly_t poly, gr_ctx_t ctx)
{
    return _gr_mat_companion(res, poly->coeffs, poly->length, ctx);
}

int
_gr_mat_companion_fraction(gr_mat_t res_num, gr_ptr res_den, gr_srcptr poly, slong len, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, n;
    slong sz = ctx->sizeof_elem;

    n = gr_mat_nrows(res_num, ctx);

    if (n == 0)
        return gr_one(res_den, ctx);

    status |= gr_set(res_den, GR_ENTRY(poly, n, sz), ctx);

    for (i = 0; i < n - 1; i++) {
        for (j = 0; j < n; j++) {
            if (i + 1 == j)
                status |= gr_set(gr_mat_entry_ptr(res_num, i, j, ctx), res_den, ctx);
            else
                status |= gr_zero(gr_mat_entry_ptr(res_num, i, j, ctx), ctx);
        }
    }

    for (j = 0; j < n; j++)
        status |= gr_neg(gr_mat_entry_ptr(res_num, n - 1, j, ctx), GR_ENTRY(poly, j, sz), ctx);

    return status;
}

int
gr_mat_companion_fraction(gr_mat_t res_num, gr_ptr res_den, const gr_poly_t poly, gr_ctx_t ctx)
{
    return _gr_mat_companion_fraction(res_num, res_den, poly->coeffs, poly->length, ctx);
}
