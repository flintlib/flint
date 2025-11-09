/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

/* To do: over rings with coefficent explosion, instead of checking that
   a dot product or matrix product is zero, compute two half-length
   products and check that they are neg-equal. */

/* Check that the square submatrix [a,b) x [a,b) in A*A^T is the identity matrix. */
static truth_t
gr_mat_is_orthogonal_recursive(const gr_mat_t A, gr_ptr tmp, slong a, slong b, int unit, gr_ctx_t ctx)
{
    slong i, j;
    truth_t ok, res = T_TRUE;
    slong ncols = A->c;
    int status;
    slong sz = ctx->sizeof_elem;

    if (b - a <= 10)
    {
        for (i = a; i < b; i++)
        {
            for (j = a; j <= (unit ? i : i - 1); j++)
            {
                status = _gr_vec_dot(tmp, NULL, 0, GR_MAT_ENTRY(A, i, 0, sz),
                                                 GR_MAT_ENTRY(A, j, 0, sz), ncols, ctx);

                ok = (status == GR_SUCCESS) ?
                    ((i == j) ? gr_is_one(tmp, ctx) : gr_is_zero(tmp, ctx)) : T_UNKNOWN;
                res = truth_and(res, ok);
                if (res == T_FALSE)
                    return res;
            }
        }
    }
    else
    {
        slong m = a + (b - a) / 2;

        res = truth_and(res, gr_mat_is_orthogonal_recursive(A, tmp, a, m, unit, ctx));
        if (res == T_FALSE)
            return res;

        res = truth_and(res, gr_mat_is_orthogonal_recursive(A, tmp, m, b, unit, ctx));
        if (res == T_FALSE)
            return res;

        /* Check that [m,b) x [a,m) in A*A^T is the zero matrix. */
        gr_mat_t A1, A2, A2t, T;

        gr_mat_window_init(A1, A, m, 0, b, ncols, ctx);
        gr_mat_window_init(A2, A, a, 0, m, ncols, ctx);

        /* There is no transposed gr_mat_mul yet, so do a shallow transpose */
        slong r = A2->r;
        slong c = A2->c;
        A2t->r = c;
        A2t->c = r;
        A2t->stride = A2t->c;
        A2t->entries = GR_TMP_ALLOC(sz * r * c);
        for (i = 0; i < r; i++)
            for (j = 0; j < c; j++)
                gr_set_shallow(GR_MAT_ENTRY(A2t, j, i, sz), GR_MAT_ENTRY(A2, i, j, sz), ctx);

        gr_mat_init(T, b - m, m - a, ctx);

        if (gr_mat_mul(T, A1, A2t, ctx) == GR_SUCCESS)
            res = gr_mat_is_zero(T, ctx);
        else
            res = T_UNKNOWN;

        GR_TMP_FREE(A2t->entries, sz * r * c);
        gr_mat_clear(T, ctx);
    }

    return res;
}

static truth_t
_gr_mat_is_row_orthogonal(const gr_mat_t A, int unit, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(A, ctx);
    truth_t res;
    gr_ptr t;

    if (n == 0)
        return T_TRUE;

    if (unit && (n > gr_mat_ncols(A, ctx) && gr_ctx_is_zero_ring(ctx) == T_FALSE))
        return T_FALSE;

    GR_TMP_INIT(t, ctx);
    res = gr_mat_is_orthogonal_recursive(A, t, 0, n, unit, ctx);
    GR_TMP_CLEAR(t, ctx);

    return res;
}

static truth_t
_gr_mat_is_col_orthogonal(const gr_mat_t A, int unit, gr_ctx_t ctx)
{
    slong n = gr_mat_ncols(A, ctx);
    truth_t res;
    gr_ptr t;
    gr_mat_t At;

    if (n == 0)
        return T_TRUE;

    if (unit && (n > gr_mat_nrows(A, ctx) && gr_ctx_is_zero_ring(ctx) == T_FALSE))
        return T_FALSE;

    /* todo: let gr_mat_is_orthogonal_recursive exploit the fact that
       we have (At^T) */
    GR_MAT_TMP_INIT_SHALLOW_TRANSPOSE(At, A, ctx);
    GR_TMP_INIT(t, ctx);
    res = gr_mat_is_orthogonal_recursive(At, t, 0, n, unit, ctx);
    GR_TMP_CLEAR(t, ctx);
    GR_MAT_TMP_CLEAR_SHALLOW_TRANSPOSE(At, ctx);

    return res;
}

truth_t
gr_mat_is_orthogonal(const gr_mat_t A, gr_ctx_t ctx)
{
    slong n = gr_mat_nrows(A, ctx);

    if (n != gr_mat_ncols(A, ctx))
        return T_FALSE;

    return _gr_mat_is_row_orthogonal(A, 1, ctx);
}

truth_t
gr_mat_is_row_orthogonal(const gr_mat_t A, gr_ctx_t ctx)
{
    return _gr_mat_is_row_orthogonal(A, 0, ctx);
}

truth_t
gr_mat_is_row_orthonormal(const gr_mat_t A, gr_ctx_t ctx)
{
    return _gr_mat_is_row_orthogonal(A, 1, ctx);
}

truth_t
gr_mat_is_col_orthogonal(const gr_mat_t A, gr_ctx_t ctx)
{
    return _gr_mat_is_col_orthogonal(A, 0, ctx);
}

truth_t
gr_mat_is_col_orthonormal(const gr_mat_t A, gr_ctx_t ctx)
{
    return _gr_mat_is_col_orthogonal(A, 1, ctx);
}

