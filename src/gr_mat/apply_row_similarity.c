/*
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

/* Add t to the entries x[0], ..., x[len - 1]. */
static int
_gr_vec_add_scalar_small(gr_ptr x, slong len, gr_srcptr t, gr_method_binary_op add, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong k, sz = ctx->sizeof_elem;

    if (len <= 2)
    {
        for (k = 0; k < len; k++)
            status |= add(GR_ENTRY(x, k, sz), GR_ENTRY(x, k, sz), t, ctx);
    }
    else
    {
        status |= _gr_vec_add_scalar(x, x, len, t, ctx);
    }

    return status;
}

/* Version for exact arithmetic: A[i][r] * d is computed once per row i, and
   (when there are more than two of them) the rows j are summed before
   multiplying by d. By distributivity this gives the same result as the
   entrywise version below. */
static int
_gr_mat_apply_row_similarity_exact(gr_mat_t A, slong r, gr_ptr d, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    int status = GR_SUCCESS;
    slong n = A->r, i, j, j0, m;
    slong sz = ctx->sizeof_elem;
    gr_ptr t, s;

    /* number of rows (or columns) j with j < r - 1 or j > r; m >= 1 */
    m = FLINT_MAX(r - 1, 0) + (n - r - 1);

    /* the first such index j */
    j0 = (r - 1 > 0) ? 0 : r + 1;

    /* columns j += d * column r */
    GR_TMP_INIT(t, ctx);

    for (i = 0; i < n; i++)
    {
        status |= mul(t, GR_MAT_ENTRY(A, i, r, sz), d, ctx);

        if (r - 1 > 0)
            status |= _gr_vec_add_scalar_small(GR_MAT_ENTRY(A, i, 0, sz), r - 1, t, add, ctx);

        if (n - r - 1 > 0)
            status |= _gr_vec_add_scalar_small(GR_MAT_ENTRY(A, i, r + 1, sz), n - r - 1, t, add, ctx);
    }

    GR_TMP_CLEAR(t, ctx);

    /* row r -= d * rows j */
    if (m <= 2)
    {
        for (j = j0; j < n; j++)
        {
            if (j == r - 1 || j == r)
                continue;

            status |= _gr_vec_submul_scalar(GR_MAT_ENTRY(A, r, 0, sz), GR_MAT_ENTRY(A, j, 0, sz), n, d, ctx);
        }
    }
    else
    {
        GR_TMP_INIT_VEC(s, n, ctx);

        status |= _gr_vec_set(s, GR_MAT_ENTRY(A, j0, 0, sz), n, ctx);

        for (j = j0 + 1; j < n; j++)
        {
            if (j == r - 1 || j == r)
                continue;

            status |= _gr_vec_add(s, s, GR_MAT_ENTRY(A, j, 0, sz), n, ctx);
        }

        status |= _gr_vec_mul_scalar(s, s, n, d, ctx);
        status |= _gr_vec_sub(GR_MAT_ENTRY(A, r, 0, sz), GR_MAT_ENTRY(A, r, 0, sz), s, n, ctx);

        GR_TMP_CLEAR_VEC(s, n, ctx);
    }

    return status;
}

int gr_mat_apply_row_similarity(gr_mat_t A, slong r, gr_ptr d, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong n = A->r, i, j;
    slong sz = ctx->sizeof_elem;

    if (A->r != A->c || r < 0 || r >= A->r)
        return GR_DOMAIN;

    /* nothing to do if there are no rows (or columns) j with j < r - 1 or j > r */
    if (r <= 1 && r == n - 1)
        return GR_SUCCESS;

    if (gr_ctx_is_exact(ctx) == T_TRUE)
        return _gr_mat_apply_row_similarity_exact(A, r, d, ctx);

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < r - 1; j++)
            status |= gr_addmul(GR_MAT_ENTRY(A, i, j, sz), GR_MAT_ENTRY(A, i, r, sz), d, ctx);

        for (j = r + 1; j < n; j++)
            status |= gr_addmul(GR_MAT_ENTRY(A, i, j, sz), GR_MAT_ENTRY(A, i, r, sz), d, ctx);
    }

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < r - 1; j++)
            status |= gr_submul(GR_MAT_ENTRY(A, r, i, sz), GR_MAT_ENTRY(A, j, i, sz), d, ctx);

        for (j = r + 1; j < n; j++)
            status |= gr_submul(GR_MAT_ENTRY(A, r, i, sz), GR_MAT_ENTRY(A, j, i, sz), d, ctx);
    }

    return status;
}
