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

int
gr_mat_rref_lu(slong * res_rank, gr_mat_t R, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j, k, n, rank;
    slong *pivots;
    slong *nonpivots;
    slong *P;
    gr_mat_t U, V;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    if (gr_mat_is_zero(A, ctx) == T_TRUE)
    {
        *res_rank = 0;
        return GR_SUCCESS;
    }

    /* Todo: fast path for nrows == 1 */

    n = A->c;

    P = _perm_init(gr_mat_nrows(A, ctx));
    status = gr_mat_lu(&rank, P, R, A, 0, ctx);
    _perm_clear(P);

    if (status != GR_SUCCESS)
        return status;

    if (rank == 0)
    {
        *res_rank = 0;
        return GR_SUCCESS;
    }

    /* Clear L */
    for (i = 0; i < A->r; i++)
        for (j = 0; j < FLINT_MIN(i, rank); j++)
            status |= gr_zero(GR_MAT_ENTRY(R, i, j, sz), ctx);

    /* We now reorder U to proper upper triangular form U | V
       with U full-rank triangular, set V = U^(-1) V, and then
       put the column back in the original order.

       An improvement for some matrices would be to compress V by
       discarding columns containing nothing but zeros. */

    gr_mat_init(U, rank, rank, ctx);
    gr_mat_init(V, rank, n - rank, ctx);

    pivots = flint_malloc(sizeof(slong) * rank);
    nonpivots = flint_malloc(sizeof(slong) * (n - rank));

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

    for (i = 0; i < rank; i++)
        for (j = 0; j <= i; j++)
            status |= gr_set(GR_MAT_ENTRY(U, j, i, sz), GR_MAT_ENTRY(R, j, pivots[i], sz), ctx);

    for (i = 0; i < n - rank; i++)
        for (j = 0; j < rank; j++)
            status |= gr_set(GR_MAT_ENTRY(V, j, i, sz), GR_MAT_ENTRY(R, j, nonpivots[i], sz), ctx);

    status |= gr_mat_nonsingular_solve_triu(V, U, V, 0, ctx);

    /* Clear pivot columns */
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
        {
            if (i == j)
                status |= gr_one(GR_MAT_ENTRY(R, j, pivots[i], sz), ctx);
            else
                status |= gr_zero(GR_MAT_ENTRY(R, j, pivots[i], sz), ctx);
        }
    }

    /* Write back the actual content */
    for (i = 0; i < n - rank; i++)
        for (j = 0; j < rank; j++)
            status |= gr_set(GR_MAT_ENTRY(R, j, nonpivots[i], sz), GR_MAT_ENTRY(V, j, i, sz), ctx);

cleanup1:
    gr_mat_clear(U, ctx);
    gr_mat_clear(V, ctx);

    flint_free(pivots);
    flint_free(nonpivots);

    *res_rank = rank;
    return status;
}
