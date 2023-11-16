/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "perm.h"
#include "ca_mat.h"

int
ca_mat_rref_lu(slong * res_rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j, k, n, rank;
    slong *pivots;
    slong *nonpivots;
    slong *P;
    ca_mat_t U, V;
    int success;

    if (ca_mat_check_is_zero(A, ctx) == T_TRUE)
    {
        *res_rank = 0;
        return 1;
    }

    /* Todo: fast path for nrows == 1 */

    n = A->c;

    P = _perm_init(ca_mat_nrows(A));
    success = ca_mat_lu(&rank, P, R, A, 0, ctx);
    _perm_clear(P);

    if (!success)
        return 0;

    if (rank == 0)
    {
        *res_rank = 0;
        return 1;
    }

    /* Clear L */
    for (i = 0; i < A->r; i++)
        for (j = 0; j < FLINT_MIN(i, rank); j++)
            ca_zero(ca_mat_entry(R, i, j), ctx);

    /* We now reorder U to proper upper triangular form U | V
       with U full-rank triangular, set V = U^(-1) V, and then
       put the column back in the original order.

       An improvement for some matrices would be to compress V by
       discarding columns containing nothing but zeros. */

    ca_mat_init(U, rank, rank, ctx);
    ca_mat_init(V, rank, n - rank, ctx);

    pivots = flint_malloc(sizeof(slong) * rank);
    nonpivots = flint_malloc(sizeof(slong) * (n - rank));

    for (i = j = k = 0; i < rank; i++)
    {
        while (1)
        {
            /* Todo: this should not be T_UNKNOWN. Should we save
               the pivot data in the lu algorithm? */
            truth_t is_zero = ca_check_is_zero(ca_mat_entry(R, i, j), ctx);

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
                success = 0;
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
            ca_set(ca_mat_entry(U, j, i), ca_mat_entry(R, j, pivots[i]), ctx);

    for (i = 0; i < n - rank; i++)
        for (j = 0; j < rank; j++)
            ca_set(ca_mat_entry(V, j, i), ca_mat_entry(R, j, nonpivots[i]), ctx);

    ca_mat_solve_triu(V, U, V, 0, ctx);

    /* Clear pivot columns */
    for (i = 0; i < rank; i++)
    {
        for (j = 0; j <= i; j++)
        {
            if (i == j)
                ca_one(ca_mat_entry(R, j, pivots[i]), ctx);
            else
                ca_zero(ca_mat_entry(R, j, pivots[i]), ctx);
        }
    }

    /* Write back the actual content */
    for (i = 0; i < n - rank; i++)
        for (j = 0; j < rank; j++)
            ca_set(ca_mat_entry(R, j, nonpivots[i]), ca_mat_entry(V, j, i), ctx);

cleanup1:
    ca_mat_clear(U, ctx);
    ca_mat_clear(V, ctx);

    flint_free(pivots);
    flint_free(nonpivots);

    *res_rank = rank;
    return 1;
}
