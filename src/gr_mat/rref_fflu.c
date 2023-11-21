/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "gr_mat.h"

static int
_gr_mat_rref_fflu(slong * res_rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, int divided, gr_ctx_t ctx)
{
    slong i, j, k, m, n, rank;
    slong *pivots;
    slong *nonpivots;
    slong *P;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (gr_mat_is_zero(A, ctx) == T_TRUE)
    {
        *res_rank = 0;
        status |= gr_one(den, ctx);
        return status;
    }

    P = _perm_init(gr_mat_nrows(A, ctx));
    status |= gr_mat_fflu(&rank, P, R, den, A, 0, ctx);
    _perm_clear(P);

    if (status != GR_SUCCESS)
        return status;

    m = gr_mat_nrows(R, ctx);
    n = gr_mat_ncols(R, ctx);

    /* clear bottom */
    for (i = rank; i < m; i++)
        for (j = 0; j < n; j++)
            status |= gr_zero(GR_MAT_ENTRY(R, i, j, sz), ctx);

    /* Convert row echelon form to reduced row echelon form */
    if (rank > 1)
    {
        gr_ptr t, u;

        GR_TMP_INIT2(t, u, ctx);

        pivots = flint_malloc(sizeof(slong) * n);
        nonpivots = pivots + rank;

        for (i = j = k = 0; i < rank; i++)
        {
            while (1)
            {
                /* Todo: this should not be T_UNKNOWN. Should we save
                   the pivot data in the lu algorithm? */
                truth_t is_zero = gr_is_zero(GR_MAT_ENTRY(R, i, j, sz), ctx);

                if (is_zero == T_FALSE)
                {
                    break;
                }
                else if (is_zero == T_TRUE)
                {
                    nonpivots[k] = j;
                    k++;
                    j++;
                }
                else
                {
                    status = GR_UNABLE;
                    goto cleanup1;
                }
            }
            pivots[i] = j;
            j++;
        }
        while (k < n - rank)
        {
            nonpivots[k] = j;
            k++;
            j++;
        }

        for (k = 0; k < n - rank; k++)
        {
            for (i = rank - 2; i >= 0; i--)
            {
                status |= gr_mul(t, den, GR_MAT_ENTRY(R, i, nonpivots[k], sz), ctx);

                for (j = i + 1; j < rank; j++)
                {
                    /* should be submul */
                    status |= gr_mul(u, GR_MAT_ENTRY(R, i, pivots[j], sz), GR_MAT_ENTRY(R, j, nonpivots[k], sz), ctx);
                    status |= gr_sub(t, t, u, ctx);
                }

                /* should be divexact */
                status |= gr_div(GR_MAT_ENTRY(R, i, nonpivots[k], sz), t, GR_MAT_ENTRY(R, i, pivots[i], sz), ctx);
            }
        }

        /* clear pivot columns */
        for (i = 0; i < rank; i++)
        {
            for (j = 0; j < rank; j++)
            {
                if (i == j)
                    if (divided)
                        status |= gr_one(GR_MAT_ENTRY(R, j, pivots[i], sz), ctx);
                    else
                        status |= gr_set(GR_MAT_ENTRY(R, j, pivots[i], sz), den, ctx);
                else
                    status |= gr_zero(GR_MAT_ENTRY(R, j, pivots[i], sz), ctx);
            }
        }

        /* divide out denominator */
        if (divided && gr_is_one(den, ctx) != T_TRUE)
            for (i = 0; i < rank; i++)
                for (j = 0; j < n - rank; j++)
                    status |= gr_div(GR_MAT_ENTRY(R, i, nonpivots[j], sz), GR_MAT_ENTRY(R, i, nonpivots[j], sz), den, ctx);

cleanup1:
        flint_free(pivots);
        GR_TMP_CLEAR2(t, u, ctx);
    }
    else if (rank == 1 && divided && gr_is_one(den, ctx) != T_TRUE)
    {
        for (i = 0; i < n; i++)
            status |= gr_div(GR_MAT_ENTRY(R, 0, i, sz), GR_MAT_ENTRY(R, 0, i, sz), den, ctx);
    }

    *res_rank = rank;

    return status;
}

int
gr_mat_rref_den_fflu(slong * res_rank, gr_mat_t R, gr_ptr den, const gr_mat_t A, gr_ctx_t ctx)
{
    return _gr_mat_rref_fflu(res_rank, R, den, A, 0, ctx);
}

int
gr_mat_rref_fflu(slong * res_rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
{
    int status;
    gr_ptr den;
    GR_TMP_INIT(den, ctx);
    status = _gr_mat_rref_fflu(res_rank, R, den, A, 1, ctx);
    GR_TMP_CLEAR(den, ctx);
    return status;
}
