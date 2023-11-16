/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

/* Defined in t-rref_den_fflu.c, t-rref_fflu.c and t-rref_lu.c */
#ifndef gr_mat_randrowops
#define gr_mat_randrowops gr_mat_randrowops
/* todo: make a function */
int
gr_mat_randrowops(gr_mat_t mat, flint_rand_t state, slong count, gr_ctx_t ctx)
{
    slong c, i, j, k;
    slong m = mat->r;
    slong n = mat->c;
    slong sz = ctx->sizeof_elem;
    int status = GR_SUCCESS;

    if (mat->r == 0 || mat->c == 0)
        return GR_SUCCESS;

    for (c = 0; c < count; c++)
    {
        if ((i = n_randint(state, m)) == (j = n_randint(state, m)))
            continue;
        if (n_randint(state, 2))
            for (k = 0; k < n; k++)
                status |= gr_add(GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, i, k, sz), ctx);
        else
            for (k = 0; k < n; k++)
                status |= gr_sub(GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, j, k, sz), GR_MAT_ENTRY(mat, i, k, sz), ctx);
    }

    return status;
}
#endif

TEST_GR_FUNCTION_START(gr_mat_rref_den_fflu, state, count_success, count_unable, count_domain)
{
    slong iter;

    /* Check that random row/column operations preserve rank */
    for (iter = 0; iter < 10000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A, B, R, R2;
        gr_ptr den;
        slong rank1, rank2, r, c;
        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        r = n_randint(state, 6);
        c = n_randint(state, 6);

        gr_mat_init(A, r, c, ctx);
        gr_mat_init(B, r, c, ctx);
        gr_mat_init(R, r, c, ctx);
        gr_mat_init(R2, r, c, ctx);
        den = gr_heap_init(ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_set(B, A, ctx);
        status |= gr_mat_randrowops(B, state, 1 + n_randint(state, 20), ctx);

        status |= gr_mat_rref_fflu(&rank1, R, A, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_mat_rref_den_fflu(&rank2, R2, den, B, ctx);
            status |= gr_mat_div_scalar(R2, R2, den, ctx);

            if (status == GR_SUCCESS)
            {
                if (rank1 != rank2 || gr_mat_equal(R, R2, ctx) == T_FALSE)
                {
                    flint_printf("FAIL:\n");
                    gr_ctx_println(ctx);
                    flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                    flint_printf("B: "); gr_mat_print(B, ctx); flint_printf("\n");
                    flint_printf("R: "); gr_mat_print(R, ctx); flint_printf("\n");
                    flint_printf("R2: "); gr_mat_print(R2, ctx); flint_printf("\n");
                    flint_printf("rank = %wd, %wd\n\n", rank1, rank2);
                    flint_abort();
                }
            }
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(R, ctx);
        gr_mat_clear(R2, ctx);
        gr_heap_clear(den, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
