/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ca_mat.h"

int
ca_mat_rref_fflu(slong * res_rank, ca_mat_t R, const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j, k, m, n, rank;
    slong *pivots, *nonpivots;
    ca_t den;
    int success;

    ca_init(den, ctx);
    success = ca_mat_fflu(&rank, NULL, R, den, A, 0, ctx);
    if (!success)
    {
        ca_clear(den, ctx);
        return 0;
    }

    m = ca_mat_nrows(R);
    n = ca_mat_ncols(R);

    /* clear bottom */
    for (i = rank; i < m; i++)
        for (j = 0; j < n; j++)
            ca_zero(ca_mat_entry(R, i, j), ctx);

    /* Convert row echelon form to reduced row echelon form */
    if (rank > 1)
    {
        ca_t t, u;

        ca_init(t, ctx);
        ca_init(u, ctx);

        pivots = flint_malloc(sizeof(slong) * n);
        nonpivots = pivots + rank;

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

        for (k = 0; k < n - rank; k++)
        {
            for (i = rank - 2; i >= 0; i--)
            {
                ca_mul(t, den, ca_mat_entry(R, i, nonpivots[k]), ctx);

                for (j = i + 1; j < rank; j++)
                {
                    ca_mul(u, ca_mat_entry(R, i, pivots[j]), ca_mat_entry(R, j, nonpivots[k]), ctx);
                    ca_sub(t, t, u, ctx);
                }

                ca_div(ca_mat_entry(R, i, nonpivots[k]), t, ca_mat_entry(R, i, pivots[i]), ctx);
            }
        }

        /* clear pivot columns */
        for (i = 0; i < rank; i++)
        {
            for (j = 0; j < rank; j++)
            {
                if (i == j)
                    ca_one(ca_mat_entry(R, j, pivots[i]), ctx);
                else
                    ca_zero(ca_mat_entry(R, j, pivots[i]), ctx);
            }
        }

        /* divide out denominator */
        if (ca_check_is_one(den, ctx) != T_TRUE)
            for (i = 0; i < rank; i++)
                for (j = 0; j < n - rank; j++)
                    ca_div(ca_mat_entry(R, i, nonpivots[j]), ca_mat_entry(R, i, nonpivots[j]), den, ctx);
cleanup1:
        flint_free(pivots);
        ca_clear(t, ctx);
        ca_clear(u, ctx);
    }
    else if (rank == 1 && ca_check_is_one(den, ctx) != T_TRUE)
    {
        for (i = 0; i < n; i++)
            ca_div(ca_mat_entry(R, 0, i), ca_mat_entry(R, 0, i), den, ctx);
    }

    ca_clear(den, ctx);
    *res_rank = rank;

    return success;
}
