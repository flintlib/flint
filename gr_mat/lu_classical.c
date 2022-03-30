/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

/*
static truth_t
ca_check_is_zero_fast(const ca_t x, ca_ctx_t ctx)
{
    if (CA_IS_QQ(x, ctx))
    {
        return fmpq_is_zero(CA_FMPQ(x)) ? T_TRUE : T_FALSE;
    }
    else
    {
        return T_UNKNOWN;
    }
}

static truth_t ca_check_is_zero_and_simplify(ca_t x, ca_ctx_t ctx)
{
    truth_t res = ca_check_is_zero_fast(x, ctx);

    if (res == T_UNKNOWN)
    {
        res = ca_check_is_zero(x, ctx);

        if (res == T_TRUE)
            ca_zero(x, ctx);
    }

    return res;
}
*/

int
gr_cmp_repr(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    return 0;
}

int
gr_mat_find_pivot(slong * pivot_row, gr_mat_t mat, slong start_row, slong end_row, slong column, gr_ctx_t ctx)
{
    slong best_row, i;
    int unknown;
    truth_t is_zero;
    slong sz;

    if (end_row <= start_row)
        flint_abort();

/* todo: for ca-type rings */
#if 0
    /* First find the simplest element that is not trivially zero. With high probability
       this will actually be nonzero. */
    best_row = -1;

    for (i = start_row; i < end_row; i++)
    {
        is_zero = ca_check_is_zero_fast(ca_mat_entry(mat, i, column), ctx);

        if (is_zero != T_TRUE)
        {
            if (best_row == -1 || ca_cmp_repr(ca_mat_entry(mat, i, column), ca_mat_entry(mat, best_row, column), ctx) < 0)
            {
                best_row = i;
            }
        }
    }

    if (best_row != -1)
    {
        is_zero = ca_check_is_zero_and_simplify(ca_mat_entry(mat, best_row, column), ctx);

        if (is_zero == T_FALSE)
        {
            *pivot_row = best_row;
            return T_TRUE;
        }
    }
#endif

    /* If the above failed, go through all elements again and do more expensive checks.
       Todo: 1) support in-place simplifications.
             2) consider sorting above and traversing all entries in order of
                simplicity (not just the simplest element). */

    best_row = -1;
    unknown = 0;

    sz = ctx->sizeof_elem;

    for (i = start_row; i < end_row; i++)
    {
        is_zero = gr_is_zero(GR_MAT_ENTRY(mat, i, column, sz), ctx);

        if (is_zero != T_UNKNOWN)
        {
            if (is_zero == T_FALSE)
            {
                if (best_row == -1 || gr_cmp_repr(GR_MAT_ENTRY(mat, i, column, sz), GR_MAT_ENTRY(mat, best_row, column, sz), ctx) < 0)
                {
                    best_row = i;
                }
            }
        }
        else
        {
            unknown = 1;
        }
    }

    if (best_row == -1)
    {
        *pivot_row = -1;
        if (unknown)
            return GR_UNABLE;
        else
            return GR_DOMAIN;
    }
    else
    {
        *pivot_row = best_row;
        return GR_SUCCESS;
    }
}

void
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
gr_mat_lu_classical(slong * res_rank, slong * P, gr_mat_t LU, const gr_mat_t A, int full_rank_check, gr_ctx_t ctx)
{
    gr_ptr d, e;
    gr_ptr * a;
    slong i, j, m, n, r, rank, row, col, sz;
    int status;
    int pivot_status;
    GR_TMP_START;

    if (gr_mat_is_empty(A, ctx) == T_TRUE)
    {
        *res_rank = 0;
        return GR_SUCCESS;
    }

    GR_TMP_INIT2(d, e, ctx);

    m = gr_mat_nrows(A, ctx);
    n = gr_mat_ncols(A, ctx);
    sz = ctx->sizeof_elem;

    gr_mat_set(LU, A, ctx);

    a = LU->rows;

    rank = row = col = 0;
    for (i = 0; i < m; i++)
        P[i] = i;

    gr_init(d, ctx);
    gr_init(e, ctx);

    status = GR_SUCCESS;

    while (row < m && col < n)
    {
        pivot_status = gr_mat_find_pivot(&r, LU, row, m, col, ctx);

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
            if (full_rank_check)
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
            status |= _gr_vec_scalar_addmul(GR_ENTRY(a[j], col + 1, sz), GR_ENTRY(a[row], col + 1, sz), n - col - 1, e, ctx);
            status |= gr_zero(GR_ENTRY(a[j], col, sz), ctx);
            status |= gr_neg(GR_ENTRY(a[j], rank - 1, sz), e, ctx);
        }

        row++;
        col++;
    }

    GR_TMP_CLEAR2(d, e, ctx);
    GR_TMP_END;

    *res_rank = rank;
    return status;
}
