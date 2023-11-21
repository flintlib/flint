/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "gr_mat.h"

/* todo: rewrite scalar operations using dots? */
/* todo: other algorithms (+ memory efficiency) */

int
_gr_mat_gr_poly_evaluate(gr_mat_t y, gr_srcptr poly,
    slong len, const gr_mat_t x, gr_ctx_t ctx)
{
    slong i, j, m, r, n;
    gr_mat_struct * xs;
    gr_mat_t s, t;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    n = gr_mat_nrows(x, ctx);

    if (n != gr_mat_ncols(x, ctx))
        return GR_DOMAIN;

    if (len == 0)
    {
        return gr_mat_zero(y, ctx);
    }

    if (len == 1)
    {
        return gr_mat_set_scalar(y, poly, ctx);
    }

    if (len == 2)
    {
        status |= gr_mat_mul_scalar(y, x, GR_ENTRY(poly, 1, sz), ctx);
        status |= gr_mat_add_scalar(y, y, poly, ctx);
        return status;
    }

    m = n_sqrt(len) + 1;
    r = (len + m - 1) / m;

    xs = flint_malloc(sizeof(gr_mat_struct) * (m + 1));
    for (i = 0; i <= m; i++)
    {
        gr_mat_init(xs + i, n, n, ctx);

        if (i == 0)
            status |= gr_mat_one(xs + i, ctx);
        else if (i == 1)
            status |= gr_mat_set(xs + i, x, ctx);
        else
            status |= gr_mat_mul(xs + i, xs + i - 1, x, ctx);
            /* todo: squaring when that is better */
    }

    gr_mat_init(s, n, n, ctx);
    gr_mat_init(t, n, n, ctx);

    status |= gr_mat_set_scalar(y, GR_ENTRY(poly, (r - 1) * m, sz), ctx);
    for (j = 1; (r - 1) * m + j < len; j++)
        status |= gr_mat_addmul_scalar(y, xs + j, GR_ENTRY(poly, (r - 1) * m + j, sz), ctx);

    for (i = r - 2; i >= 0; i--)
    {
        status |= gr_mat_set_scalar(s, GR_ENTRY(poly, i * m, sz), ctx);
        for (j = 1; j < m; j++)
            status |= gr_mat_addmul_scalar(s, xs + j, GR_ENTRY(poly, i * m + j, sz), ctx);

        status |= gr_mat_mul(y, y, xs + m, ctx);
        status |= gr_mat_add(y, y, s, ctx);
    }

    for (i = 0; i <= m; i++)
        gr_mat_clear(xs + i, ctx);
    flint_free(xs);
    gr_mat_clear(s, ctx);
    gr_mat_clear(t, ctx);

    return status;
}

int
gr_mat_gr_poly_evaluate(gr_mat_t res, const gr_poly_t f, const gr_mat_t a, gr_ctx_t ctx)
{
    return _gr_mat_gr_poly_evaluate(res, f->coeffs, f->length, a, ctx);
}
