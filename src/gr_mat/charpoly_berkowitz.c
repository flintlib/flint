/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

int
_gr_mat_charpoly_berkowitz(gr_ptr cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    slong n = mat->r;
    int status;

    status = GR_SUCCESS;

    if (n == 0)
    {
        status |= gr_one(cp, ctx);
    }
    else if (n == 1)
    {
        status |= gr_neg(cp, mat->rows[0], ctx);
        status |= gr_one(GR_ENTRY(cp, 1, sz), ctx);
    }
    else if (n == 2)
    {
        status |= gr_mat_det_cofactor(cp, mat, ctx);
        status |= gr_add(GR_ENTRY(cp, 1, sz), GR_MAT_ENTRY(mat, 0, 0, sz), GR_MAT_ENTRY(mat, 1, 1, sz), ctx);
        status |= gr_neg(GR_ENTRY(cp, 1, sz), GR_ENTRY(cp, 1, sz), ctx);
        status |= gr_one(GR_ENTRY(cp, 2, sz), ctx);
    }
    else
    {
        slong i, k, t;
        gr_ptr a, A, s;

        GR_TMP_INIT_VEC(a, n * n, ctx);
        A = GR_ENTRY(a, (n - 1) * n, sz);

        status |= _gr_vec_zero(cp, n + 1, ctx);
        status |= gr_neg(cp, GR_MAT_ENTRY(mat, 0, 0, sz), ctx);

        for (t = 1; t < n; t++)
        {
            for (i = 0; i <= t; i++)
            {
                status |= gr_set(GR_ENTRY(a, 0 * n + i, sz), GR_MAT_ENTRY(mat, i, t, sz), ctx);
            }

            status |= gr_set(A, GR_MAT_ENTRY(mat, t, t, sz), ctx);

            for (k = 1; k < t; k++)
            {
                for (i = 0; i <= t; i++)
                {
                    s = GR_ENTRY(a, k * n + i, sz);
                    status |= _gr_vec_dot(s, NULL, 0, mat->rows[i], GR_ENTRY(a, (k - 1) * n, sz), t + 1, ctx);
                }

                status |= gr_set(GR_ENTRY(A, k, sz), GR_ENTRY(a, k * n + t, sz), ctx);
            }

            status |= _gr_vec_dot(GR_ENTRY(A, t, sz), NULL, 0, mat->rows[t], GR_ENTRY(a, (t - 1) * n, sz), t + 1, ctx);

            for (k = 0; k <= t; k++)
            {
                status |= _gr_vec_dot_rev(GR_ENTRY(cp, k, sz), GR_ENTRY(cp, k, sz), 1, A, cp, k, ctx);
                status |= gr_sub(GR_ENTRY(cp, k, sz), GR_ENTRY(cp, k, sz), GR_ENTRY(A, k, sz), ctx);
            }
        }

        /* Shift all coefficients up by one */
        for (i = n; i > 0; i--)
            gr_swap(GR_ENTRY(cp, i, sz), GR_ENTRY(cp, (i - 1), sz), ctx);

        status |= gr_one(cp, ctx);
        status |= _gr_poly_reverse(cp, cp, n + 1, n + 1, ctx);

        GR_TMP_CLEAR_VEC(a, n * n, ctx);
    }

    return status;
}

int gr_mat_charpoly_berkowitz(gr_poly_t cp, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;

    if (mat->r != mat->c)
        return GR_DOMAIN;

    gr_poly_fit_length(cp, mat->r + 1, ctx);
    _gr_poly_set_length(cp, mat->r + 1, ctx);
    status = _gr_mat_charpoly_berkowitz(cp->coeffs, mat, ctx);
    _gr_poly_normalise(cp, ctx);   /* only needed for the zero ring */
    return status;
}
