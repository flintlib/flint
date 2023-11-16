/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gr_mat.h"

int gr_generic_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx);

/* swap (a, b) and (c, d) using (t, u) as tmp space */
static void
gr_swap2(gr_ptr a, gr_ptr b, gr_ptr c, gr_ptr d, gr_ptr t, gr_ptr u, gr_ctx_t ctx)
{
    gr_swap(t, a, ctx);
    gr_swap(u, b, ctx);
    gr_swap(a, c, ctx);
    gr_swap(b, d, ctx);
    gr_swap(c, t, ctx);
    gr_swap(d, u, ctx);
}

int
gr_mat_hessenberg_gauss(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    gr_method_binary_op mul = GR_BINARY_OP(ctx, MUL);
    gr_method_binary_op add = GR_BINARY_OP(ctx, ADD);
    gr_method_binary_op addmul = GR_BINARY_OP(ctx, ADDMUL);
    slong n, m, i, j;
    gr_ptr h, u, t;
    truth_t is_zero;
    int have_addmul;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;

    n = mat->r;
    if (n != mat->c)
        return GR_DOMAIN;

    status |= gr_mat_set(res, mat, ctx);

#define MAT_ENTRY(i, j) GR_MAT_ENTRY(res, i, j, sz)

    have_addmul = (addmul != (gr_method_binary_op) gr_generic_addmul);

#define ADDMUL(z, x, y) \
    if (have_addmul) \
    { \
        status |= addmul(z, x, y, ctx); \
    } \
    else \
    { \
        status |= mul(t, x, y, ctx); \
        status |= add(z, z, t, ctx); \
    } \

    if (n <= 2)
        return GR_SUCCESS;

    GR_TMP_INIT3(h, u, t, ctx);

    for (m = 2; m < n; m++)
    {
        i = m + 1;
        while (i <= n)
        {
            is_zero = gr_is_zero(MAT_ENTRY(i - 1, m - 1 - 1), ctx);

            if (is_zero == T_UNKNOWN)
            {
                status = GR_UNABLE;
                goto cleanup;
            }

            if (is_zero == T_FALSE)
                break;

            i += 1;
        }

        if (i != n + 1)
        {
            is_zero = gr_is_zero(MAT_ENTRY(m - 1, m - 1 - 1), ctx);

            if (is_zero == T_UNKNOWN)
            {
                status = GR_UNABLE;
                goto cleanup;
            }

            if (is_zero == T_FALSE)
                i = m;

            status |= gr_inv(h, MAT_ENTRY(i - 1, m - 1 - 1), ctx);
            if (status != GR_SUCCESS)
                goto cleanup;

            status |= gr_neg(h, h, ctx);

            if (i > m)
            {
                for (j = m - 1; j < n + 1; j++)
                    gr_swap2(MAT_ENTRY(i - 1, j - 1), MAT_ENTRY(m - 1, j - 1),
                        MAT_ENTRY(m - 1, j - 1), MAT_ENTRY(i - 1, j - 1), t, u, ctx);

                for (j = 1; j < n + 1; j++)
                    gr_swap2(MAT_ENTRY(j - 1, i - 1), MAT_ENTRY(j - 1, m - 1),
                        MAT_ENTRY(j - 1, m - 1), MAT_ENTRY(j - 1, i - 1), t, u, ctx);
            }

            for (i = m + 1; i < n + 1; i++)
            {
                is_zero = gr_is_zero(MAT_ENTRY(i - 1, m - 1 - 1), ctx);

                if (is_zero == T_UNKNOWN)
                {
                    status = GR_UNABLE;
                    goto cleanup;
                }

                if (is_zero == T_FALSE)
                {
                    status |= gr_mul(u, MAT_ENTRY(i - 1, m - 1 - 1), h, ctx);

                    for (j = m; j < n + 1; j++)
                        ADDMUL(MAT_ENTRY(i - 1, j - 1), u, MAT_ENTRY(m - 1, j - 1));

                    status |= gr_neg(u, u, ctx);

                    for (j = 1; j < n + 1; j++)
                        ADDMUL(MAT_ENTRY(j - 1, m - 1), u, MAT_ENTRY(j - 1, i - 1));

                    status |= gr_zero(MAT_ENTRY(i - 1, m - 1 - 1), ctx);
                }
            }
        }
    }

cleanup:
    GR_TMP_CLEAR3(h, u, t, ctx);

    return status;
}
