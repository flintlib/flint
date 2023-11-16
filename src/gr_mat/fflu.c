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
gr_mat_fflu(slong * res_rank, slong * P, gr_mat_t LU, gr_ptr den, const gr_mat_t A, int rank_check, gr_ctx_t ctx)
{
    gr_ptr d, e;
    gr_ptr * a;
    slong i, j, k, m, n, r, rank, row, col, sz;
    int status = GR_SUCCESS;
    int pivot_status;

    if (gr_mat_is_empty(A, ctx) == T_TRUE)
    {
        *res_rank = 0;
        return gr_one(den, ctx);
    }

    if (gr_ctx_is_integral_domain(ctx) != T_TRUE)
        return GR_UNABLE;

    GR_TMP_INIT2(d, e, ctx);

    m = gr_mat_nrows(A, ctx);
    n = gr_mat_ncols(A, ctx);
    sz = ctx->sizeof_elem;

    status |= gr_mat_set(LU, A, ctx);

    a = LU->rows;

#define ENTRY(i, j) GR_ENTRY(a[i], j, sz)

    rank = row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    gr_init(d, ctx);
    gr_init(e, ctx);
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

        /* Todo: when the pivot element is invertible, make
           use of that property.
        status |= gr_inv(d, GR_ENTRY(a[row], col, sz), ctx);
        if (status != GR_SUCCESS)
            break;
        */

        for (j = row + 1; j < m; j++)
        {
/*
            status |= gr_mul(e, ENTRY(j, col), d, ctx);
            status |= gr_neg(e, e, ctx);
            status |= _gr_vec_addmul_scalar(ENTRY(j, col + 1), ENTRY(row, col + 1), n - col - 1, e, ctx);
            status |= gr_zero(ENTRY(j, col), ctx);
            status |= gr_neg(ENTRY(j, rank - 1), e, ctx);
*/


            for (k = col + 1; k < n; k++)
            {
                status |= gr_mul(ENTRY(j, k), ENTRY(j, k), ENTRY(row, col), ctx);
                status |= gr_mul(e, ENTRY(j, col), ENTRY(row, k), ctx);
                status |= gr_sub(ENTRY(j, k), ENTRY(j, k), e, ctx);

                if (row > 0)
                {
                    status |= gr_divexact(ENTRY(j, k), ENTRY(j, k), den, ctx);

                    if (status != GR_SUCCESS)
                        goto cleanup;
                }
            }

        }

        status |= gr_set(den, ENTRY(row, col), ctx);
        row++;
        col++;
    }

cleanup:
    GR_TMP_CLEAR2(d, e, ctx);

    *res_rank = rank;
    return status;
}
