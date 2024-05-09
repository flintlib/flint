/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_poly.h"

/* todo: when adj != A, could save memory by using as temporary */
int
_gr_mat_charpoly_faddeev(gr_ptr c, gr_mat_t adj, const gr_mat_t A, gr_ctx_t ctx)
{
    gr_mat_t B, C;
    slong n, k;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    n = A->r;

    if (n == 0)
        return gr_one(c, ctx);

    if (n == 1)
    {
        status |= gr_neg(c, GR_MAT_ENTRY(A, 0, 0, sz), ctx);
        status |= gr_one(GR_ENTRY(c, 1, sz), ctx);

        if (adj != NULL)
            status |= gr_mat_one(adj, ctx);

        return status;
    }

    gr_mat_init(B, n, n, ctx);
    gr_mat_init(C, n, n, ctx);

    status |= gr_one(GR_ENTRY(c, n, sz), ctx);
    status |= gr_mat_trace(GR_ENTRY(c, n - 1, sz), A, ctx);
    status |= gr_neg(GR_ENTRY(c, n - 1, sz), GR_ENTRY(c, n - 1, sz), ctx);
    status |= gr_mat_add_scalar(B, A, GR_ENTRY(c, n - 1, sz), ctx);

    for (k = 2; k < n; k++)
    {
        status |= gr_mat_mul(C, A, B, ctx);
        status |= gr_mat_trace(GR_ENTRY(c, n - k, sz), C, ctx);
        status |= gr_div_si(GR_ENTRY(c, n - k, sz), GR_ENTRY(c, n - k, sz), -k, ctx);
        if (status != GR_SUCCESS)
            goto cleanup;

        status |= gr_mat_add_scalar(B, C, GR_ENTRY(c, n - k, sz), ctx);
    }

    status |= gr_mat_mul(C, A, B, ctx);
    status |= gr_mat_trace(c, C, ctx);
    status |= gr_div_si(c, c, -n, ctx);

    if (adj != NULL)
    {
        if (n % 2)
            status |= gr_mat_set(adj, B, ctx);
        else
            status |= gr_mat_neg(adj, B, ctx);
    }

cleanup:
    gr_mat_clear(B, ctx);
    gr_mat_clear(C, ctx);

    return status;
}

int
gr_mat_charpoly_faddeev(gr_poly_t cp, gr_mat_t adjugate, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;

    if (mat->r != mat->c)
        return GR_DOMAIN;

    gr_poly_fit_length(cp, mat->r + 1, ctx);
    _gr_poly_set_length(cp, mat->r + 1, ctx);
    status = _gr_mat_charpoly_faddeev(cp->coeffs, adjugate, mat, ctx);
    _gr_poly_normalise(cp, ctx);   /* only needed for the zero ring */
    return status;
}
