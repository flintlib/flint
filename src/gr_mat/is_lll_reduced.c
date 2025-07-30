/*
    Copyright (C) 2014 Abhinav Baid
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr.h"
#include "gr_vec.h"
#include "gr_mat.h"

truth_t
gr_mat_is_row_lll_reduced_with_removal_naive(const gr_mat_t A, gr_srcptr delta, gr_srcptr eta, gr_srcptr gs_B, slong newd, gr_ctx_t ctx)
{
    truth_t res = T_TRUE;
    int status = GR_SUCCESS;
    slong i, j, d = A->r, n = A->c;
    gr_mat_t B, mu;
    gr_ptr tmp;
    slong sz = ctx->sizeof_elem;
    int with_removal = (gs_B != NULL);

    if (d == 0 || d == 1)
        return GR_SUCCESS;

    gr_mat_init(B, d, n, ctx);
    gr_mat_init(mu, d, d, ctx);
    GR_TMP_INIT(tmp, ctx);

    status |= _gr_vec_set(GR_MAT_ENTRY(B, 0, 0, sz), GR_MAT_ENTRY(A, 0, 0, sz), n, ctx);
    /* diagonal of mu stores the squared GS norms */
    status |= _gr_vec_dot(GR_MAT_ENTRY(mu, 0, 0, sz), NULL, 0,
        GR_MAT_ENTRY(B, 0, 0, sz), GR_MAT_ENTRY(B, 0, 0, sz), n, ctx);

    if (with_removal && newd == 0)
    {
        res = gr_ge(GR_MAT_ENTRY(mu, 0, 0, sz), gs_B, ctx);
        if (res != GR_SUCCESS)
            goto cleanup;
    }

    for (i = 1; i < d; i++)
    {
        status |= _gr_vec_set(GR_MAT_ENTRY(B, i, 0, sz), GR_MAT_ENTRY(A, i, 0, sz), n, ctx);

        for (j = 0; j < i; j++)
        {
            status |= _gr_vec_dot(tmp, NULL, 0, GR_MAT_ENTRY(A, i, 0, sz), GR_MAT_ENTRY(B, j, 0, sz), n, ctx);
            status |= gr_div(GR_MAT_ENTRY(mu, i, j, sz), tmp, GR_MAT_ENTRY(mu, j, j, sz), ctx);
            if (status != GR_SUCCESS)
                goto cleanup;

            status |= _gr_vec_submul_scalar(GR_MAT_ENTRY(B, i, 0, sz), GR_MAT_ENTRY(B, j, 0, sz), n, GR_MAT_ENTRY(mu, i, j, sz), ctx);

            if (!with_removal || i < newd)
            {
                /* Check size reduction. */
                status |= gr_abs(tmp, GR_MAT_ENTRY(mu, i, j, sz), ctx);
                res = gr_le(tmp, eta, ctx);
                if (res != T_TRUE)
                    goto cleanup;
            }
        }

        status |= _gr_vec_dot(GR_MAT_ENTRY(mu, i, i, sz), NULL, 0, GR_MAT_ENTRY(B, i, 0, sz), GR_MAT_ENTRY(B, i, 0, sz), n, ctx);
        if (with_removal && i >= newd)
        {
            res = gr_ge(GR_MAT_ENTRY(mu, i, i, sz), gs_B, ctx);
            if (res != T_TRUE)
                goto cleanup;
        }

        /* Check Lovasz condition. */
        if (!with_removal || i < newd)
        {
            status |= gr_sqr(tmp, GR_MAT_ENTRY(mu, i, i - 1, sz), ctx);
            status |= gr_sub(tmp, delta, tmp, ctx);
            status |= gr_mul(tmp, tmp, GR_MAT_ENTRY(mu, i - 1, i - 1, sz), ctx);
            res = gr_le(tmp, GR_MAT_ENTRY(mu, i, i, sz), ctx);
            if (res != T_TRUE)
                goto cleanup;
        }
    }

cleanup:

    if (status != GR_SUCCESS)
        res = T_UNKNOWN;

    gr_mat_clear(B, ctx);
    gr_mat_clear(mu, ctx);
    GR_TMP_CLEAR(tmp, ctx);

    return res;
}

truth_t
gr_mat_is_row_lll_reduced_naive(const gr_mat_t A, gr_srcptr delta, gr_srcptr eta, gr_ctx_t ctx)
{
    return gr_mat_is_row_lll_reduced_with_removal_naive(A, delta, eta, NULL, 0, ctx);
}

