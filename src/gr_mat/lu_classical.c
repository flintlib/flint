/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_vec.h"
#include "gr_mat.h"

static void
_gr_mat_swap_rows(gr_mat_t mat, slong * perm, slong r, slong s, gr_ctx_t ctx)
{
    if (r != s)
    {
        gr_ptr u;
        slong t;

        if (perm != NULL)
        {
            t = perm[s];
            perm[s] = perm[r];
            perm[r] = t;
        }

        u = mat->rows[s];
        mat->rows[s] = mat->rows[r];
        mat->rows[r] = u;
    }
}

int
gr_mat_lu_classical(slong * res_rank, slong * P, gr_mat_t LU, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    gr_ptr d, e;
    gr_ptr * a;
    slong i, j, m, n, r, rank, row, col, sz;
    int status = GR_SUCCESS;
    int pivot_status;

    if (gr_mat_is_empty(A, ctx) == T_TRUE)
    {
        *res_rank = 0;
        return GR_SUCCESS;
    }

    GR_TMP_INIT2(d, e, ctx);

    m = gr_mat_nrows(A, ctx);
    n = gr_mat_ncols(A, ctx);
    sz = ctx->sizeof_elem;

    status |= gr_mat_set(LU, A, ctx);

    a = LU->rows;

    rank = row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    while (row < m && col < n)
    {
        pivot_status = gr_mat_find_nonzero_pivot(&r, LU, row, m, col, ctx);

        /* We don't know whether there is a nonzero pivot element,
           so we can't determine the rank. */
        if (pivot_status == GR_UNABLE)
        {
            status = GR_UNABLE;
            break;
        }

        /* There is certainly no nonzero pivot element. */
        if (pivot_status == GR_DOMAIN)
        {
            /* We proved that the matrix is rank-deficient,
               accomplishing the goal. */
            if (rank_check)
            {
                status = GR_SUCCESS;
                rank = 0;
                break;
            }

            /* Continue with next column. */
            col++;
            continue;
        }

        rank++;

        if (r != row)
            _gr_mat_swap_rows(LU, P, row, r, ctx);

        /* Must be able to invert pivot element. */
        status |= gr_inv(d, GR_ENTRY(a[row], col, sz), ctx);
        if (status != GR_SUCCESS)
            break;

        for (j = row + 1; j < m; j++)
        {
            status |= gr_mul(e, GR_ENTRY(a[j], col, sz), d, ctx);
            status |= gr_neg(e, e, ctx);

            if (n - col - 1 > 0)
                status |= _gr_vec_addmul_scalar(GR_ENTRY(a[j], col + 1, sz), GR_ENTRY(a[row], col + 1, sz), n - col - 1, e, ctx);

            status |= gr_zero(GR_ENTRY(a[j], col, sz), ctx);
            status |= gr_neg(GR_ENTRY(a[j], rank - 1, sz), e, ctx);
        }

        row++;
        col++;
    }

    GR_TMP_CLEAR2(d, e, ctx);

    *res_rank = rank;
    return status;
}
