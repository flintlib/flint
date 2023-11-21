/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/* todo: check that dimensions are compatible */
int
gr_mat_nonsingular_solve_fflu_precomp(gr_mat_t X, const slong * perm,
    const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;
    slong i, j, k, c, n, m;
    gr_ptr t;

    n = gr_mat_nrows(X, ctx);
    m = gr_mat_ncols(X, ctx);

    if (X == B)
    {
        gr_method_void_unary_op set_shallow = GR_VOID_UNARY_OP(ctx, SET_SHALLOW);
        /* todo: don't use malloc */
        gr_ptr tmp = flint_malloc(sz * n);

        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
                set_shallow(GR_ENTRY(tmp, i, sz), GR_MAT_ENTRY(B, perm[i], c, sz), ctx);
            for (i = 0; i < n; i++)
                set_shallow(GR_MAT_ENTRY(X, i, c, sz), GR_ENTRY(tmp, i, sz), ctx);
        }

        flint_free(tmp);
    }
    else
    {
        for (c = 0; c < m; c++)
        {
            for (i = 0; i < n; i++)
            {
                status |= gr_set(GR_MAT_ENTRY(X, i, c, sz),
                    GR_MAT_ENTRY(B, perm[i], c, sz), ctx);
            }
        }
    }

    {
        GR_TMP_INIT(t, ctx);

        /* todo: use submul */
        /* todo: use divexact? */

        for (k = 0; k < m; k++)
        {
            /* Fraction-free forward substitution */
            for (i = 0; i < n - 1; i++)
            {
                for (j = i + 1; j < n; j++)
                {
                    status |= gr_mul(GR_MAT_ENTRY(X, j, k, sz), GR_MAT_ENTRY(X, j, k, sz), GR_MAT_ENTRY(A, i, i, sz), ctx);
                    status |= gr_mul(t, GR_MAT_ENTRY(A, j, i, sz), GR_MAT_ENTRY(X, i, k, sz), ctx);
                    status |= gr_sub(GR_MAT_ENTRY(X, j, k, sz), GR_MAT_ENTRY(X, j, k, sz), t, ctx);
                    if (i > 0)
                        status |= gr_div(GR_MAT_ENTRY(X, j, k, sz), GR_MAT_ENTRY(X, j, k, sz), GR_MAT_ENTRY(A, i-1, i-1, sz), ctx);
                }
            }

            /* Fraction-free back substitution */
            for (i = n - 2; i >= 0; i--)
            {
                status |= gr_mul(GR_MAT_ENTRY(X, i, k, sz), GR_MAT_ENTRY(X, i, k, sz), GR_MAT_ENTRY(A, n-1, n-1, sz), ctx);
                for (j = i + 1; j < n; j++)
                {
                    status |= gr_mul(t, GR_MAT_ENTRY(X, j, k, sz), GR_MAT_ENTRY(A, i, j, sz), ctx);
                    status |= gr_sub(GR_MAT_ENTRY(X, i, k, sz), GR_MAT_ENTRY(X, i, k, sz), t, ctx);
                }
                status |= gr_div(GR_MAT_ENTRY(X, i, k, sz), GR_MAT_ENTRY(X, i, k, sz), GR_MAT_ENTRY(A, i, i, sz), ctx);
            }
        }

        /* status |= gr_mat_div_scalar(X, X, den, ctx); */

        GR_TMP_CLEAR(t, ctx);
    }

    return status;
}
