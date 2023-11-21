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
ca_mat_right_kernel(ca_mat_t X, const ca_mat_t A, ca_ctx_t ctx)
{
    slong i, j, k, m, n, rank, nullity;
    slong *p;
    slong *pivots;
    slong *nonpivots;
    ca_mat_t tmp;
    int success;

    m = ca_mat_nrows(A);
    n = ca_mat_ncols(A);

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));

    ca_mat_init(tmp, m, n, ctx);
    success = ca_mat_rref(&rank, tmp, A, ctx);
    nullity = n - rank;

    if (!success)
        goto cleanup;

    ca_mat_clear(X, ctx);
    ca_mat_init(X, n, nullity, ctx);

    success = 1;

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            ca_one(ca_mat_entry(X, i, i), ctx);
    }
    else if (nullity)
    {
        pivots = p;             /* length = rank */
        nonpivots = p + rank;   /* length = nullity */

        for (i = j = k = 0; i < rank; i++)
        {
            while (1)
            {
                /* Todo: this should not be T_UNKNOWN. Should we save
                   the pivot data in the lu algorithm? */
                truth_t is_zero = ca_check_is_zero(ca_mat_entry(tmp, i, j), ctx);

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
                    goto cleanup;
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

        for (i = 0; i < nullity; i++)
        {
            for (j = 0; j < rank; j++)
                ca_neg(ca_mat_entry(X, pivots[j], i), ca_mat_entry(tmp, j, nonpivots[i]), ctx);

            ca_one(ca_mat_entry(X, nonpivots[i], i), ctx);
        }
    }

cleanup:
    flint_free(p);
    ca_mat_clear(tmp, ctx);

    return success;
}
