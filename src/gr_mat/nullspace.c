/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int
gr_mat_nullspace(gr_mat_t X, const gr_mat_t A, gr_ctx_t ctx)
{
    slong i, j, k, m, n, rank, nullity;
    slong *p;
    slong *pivots;
    slong *nonpivots;
    gr_mat_t tmp;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    gr_ptr den;
    int with_den;

    m = gr_mat_nrows(A, ctx);
    n = gr_mat_ncols(A, ctx);

    p = flint_malloc(sizeof(slong) * FLINT_MAX(m, n));
    gr_mat_init(tmp, m, n, ctx);

    GR_TMP_INIT(den, ctx);

    with_den = gr_ctx_is_field(ctx) == T_FALSE && gr_ctx_is_integral_domain(ctx) == T_TRUE;

    if (with_den)
    {
        status |= gr_mat_rref_den(&rank, tmp, den, A, ctx);
    }
    else
    {
        status |= gr_mat_rref(&rank, tmp, A, ctx);
    }

    nullity = n - rank;

    if (status != GR_SUCCESS)
        goto cleanup;

    /* todo: have a resize function */
    gr_mat_clear(X, ctx);
    gr_mat_init(X, n, nullity, ctx);

    if (rank == 0)
    {
        for (i = 0; i < nullity; i++)
            status |= gr_one(GR_MAT_ENTRY(X, i, i, sz), ctx);
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
                truth_t is_zero = gr_is_zero(GR_MAT_ENTRY(tmp, i, j, sz), ctx);

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

        if (with_den)
        {
            /* if we did not keep den, equivalently den = tmp[0, pivots[0]] here */

            for (i = 0; i < nullity; i++)
            {
                for (j = 0; j < rank; j++)
                    status |= gr_neg(GR_MAT_ENTRY(X, pivots[j], i, sz), GR_MAT_ENTRY(tmp, j, nonpivots[i], sz), ctx);

                status |= gr_set(GR_MAT_ENTRY(X, nonpivots[i], + i, sz), den, ctx);
            }
        }
        else
        {
            for (i = 0; i < nullity; i++)
            {
                for (j = 0; j < rank; j++)
                    status |= gr_neg(GR_MAT_ENTRY(X, pivots[j], i, sz), GR_MAT_ENTRY(tmp, j, nonpivots[i], sz), ctx);

                status |= gr_one(GR_MAT_ENTRY(X, nonpivots[i], i, sz), ctx);
            }
        }
    }

cleanup:
    flint_free(p);
    gr_mat_clear(tmp, ctx);

    GR_TMP_CLEAR(den, ctx);

    return status;
}
