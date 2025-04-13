/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"
#include "gr_vec.h"
#include "gr_mat.h"

int
_gr_mat_companion(gr_mat_t res, gr_srcptr poly, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n = gr_mat_nrows(res, ctx);
    slong sz = ctx->sizeof_elem;
    gr_ptr c;
    gr_mat_t W;

    if (n == 0)
        return status;

    gr_mat_window_init(W, res, 0, 0, n - 1, 1, ctx);
    status |= gr_mat_zero(W, ctx);
    gr_mat_window_init(W, res, 0, 1, n - 1, n, ctx);
    status |= gr_mat_one(W, ctx);

    GR_TMP_INIT(c, ctx);
    status |= gr_inv(c, GR_ENTRY(poly, n, sz), ctx);
    status |= gr_neg(c, c, ctx);
    status |= _gr_vec_mul_scalar(gr_mat_entry_ptr(res, n - 1, 0, ctx), poly, n, c, ctx);
    GR_TMP_CLEAR(c, ctx);

    return status;
}

int
gr_mat_companion(gr_mat_t res, const gr_poly_t poly, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(res, ctx);

    if (n != gr_poly_length(poly, ctx) - 1 || n != gr_mat_ncols(res, ctx))
        return GR_DOMAIN;

    return _gr_mat_companion(res, poly->coeffs, ctx);
}

int
_gr_mat_companion_fraction(gr_mat_t res_num, gr_ptr res_den, gr_srcptr poly, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n = gr_mat_nrows(res_num, ctx);
    slong sz = ctx->sizeof_elem;
    gr_mat_t W;

    if (n == 0)
        return gr_one(res_den, ctx);

    status |= gr_set(res_den, GR_ENTRY(poly, n, sz), ctx);

    gr_mat_window_init(W, res_num, 0, 0, n - 1, 1, ctx);
    status |= gr_mat_zero(W, ctx);
    gr_mat_window_init(W, res_num, 0, 1, n - 1, n, ctx);
    status |= gr_mat_set_scalar(W, res_den, ctx);

    status |= _gr_vec_neg(gr_mat_entry_ptr(res_num, n - 1, 0, ctx), poly, n, ctx);

    return status;
}

int
gr_mat_companion_fraction(gr_mat_t res_num, gr_ptr res_den, const gr_poly_t poly, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(res_num, ctx);

    if (n != gr_poly_length(poly, ctx) - 1 || n != gr_mat_ncols(res_num, ctx))
        return GR_DOMAIN;

    return _gr_mat_companion_fraction(res_num, res_den, poly->coeffs, ctx);
}
