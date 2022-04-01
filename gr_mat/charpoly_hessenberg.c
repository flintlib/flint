/*
    Copyright (C) 2015 William Hart
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"
#include "gr_poly.h"

/* todo: optimize; rewrite using underscore methods */
int
_gr_mat_charpoly_hessenberg(gr_ptr res, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong n;
    gr_poly_struct * P;
    gr_poly_t x, v;
    gr_ptr t, u;
    slong i, m;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    GR_TMP_START;

    n = mat->r;

    gr_poly_init(x, ctx);
    gr_poly_init(v, ctx);
    GR_TMP_INIT2(t, u, ctx);

    P = flint_malloc(sizeof(gr_poly_struct) * (n + 1));
    for (i = 0; i <= n; i++)
        gr_poly_init(P + i, ctx);

    status |= gr_poly_one(P, ctx);
    status |= gr_poly_set_coeff_si(x, 1, 1, ctx);

    for (m = 1; m < n + 1; m++)
    {
        status |= gr_poly_set_scalar(v, GR_MAT_ENTRY(mat, m - 1, m - 1, sz), ctx);
        status |= gr_poly_sub(v, x, v, ctx);
        status |= gr_poly_mul(P + m + 1 - 1, v, P + m - 1, ctx);
        status |= gr_one(t, ctx);

        for (i = 1; i < m; i++)
        {
            status |= gr_mul(t, t, GR_MAT_ENTRY(mat, m - i + 1 - 1, m - i - 1, sz), ctx);
            status |= gr_poly_mul_scalar(v, P + m - i - 1, GR_MAT_ENTRY(mat, m - i - 1, m - 1, sz), ctx);
            status |= gr_poly_mul_scalar(v, v, t, ctx);
            status |= gr_poly_sub(P + m + 1 - 1, P + m + 1 - 1, v, ctx);
        }
    }

    /* for the zero ring (?) */
    status |= _gr_vec_zero(res, n + 1, ctx);
    status |= _gr_vec_set(res, (P + n + 1 - 1)->coeffs, FLINT_MIN(n + 1, (P + n + 1 - 1)->length), ctx);

    for (i = 0; i <= n; i++)
        gr_poly_clear(P + i, ctx);
    flint_free(P);

    gr_poly_clear(x, ctx);
    gr_poly_clear(v, ctx);

    GR_TMP_CLEAR2(t, u, ctx);
    GR_TMP_END;

    return status;
}

int
gr_mat_charpoly_hessenberg(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;

    if (mat->r != mat->c)
        return GR_DOMAIN;

    gr_poly_fit_length(cp, mat->r + 1, ctx);
    _gr_poly_set_length(cp, mat->r + 1, ctx);
    status = _gr_mat_charpoly_hessenberg(cp->coeffs, mat, ctx);
    _gr_poly_normalise(cp, ctx);   /* only needed for the zero ring */
    return status;
}
