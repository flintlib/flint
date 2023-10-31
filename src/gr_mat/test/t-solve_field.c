/*
    Copyright (C) 2010, 2011 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2018 Tommy Hofmann
    Copyright (C) 2020 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"

#include "fq_nmod_mat.h"

TEST_FUNCTION_START(gr_mat_solve_field, state)
{
    gr_ctx_t ctx;
    gr_mat_t A, X, X2, B, AX;
    slong i, k, m, n;
    int status;

    /* test random systems */
    for (i = 0; i < 1000; i++)
    {
        status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

/*        gr_ctx_init_fmpq(ctx);
        fmpz p = 17;
        gr_ctx_init_fq_nmod(ctx, &p, 1, "a");
*/

        m = n_randint(state, 6);
        n = n_randint(state, 6);
        k = n_randint(state, 6);

        gr_mat_init(A, m, k, ctx);
        gr_mat_init(B, m, n, ctx);
        gr_mat_init(X, k, n, ctx);
        gr_mat_init(AX, m, n, ctx);

        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_randtest(B, state, ctx);

/*
        gr_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1), ctx);
        if (n_randint(state, 2))
            gr_mat_randops(A, 1+n_randint(state, 1+m*m), state, ctx);
*/

        status |= gr_mat_solve_field(X, A, B, ctx);
        status |= gr_mat_mul(AX, A, X, ctx);

        if ((status == GR_SUCCESS && gr_mat_equal(AX, B, ctx) == T_FALSE))
        {
            flint_printf("FAIL:\n");
            flint_printf("AX != B!\n");
            flint_printf("A:\n");
            gr_mat_print(A, ctx);
            flint_printf("B:\n");
            gr_mat_print(B, ctx);
            flint_printf("X:\n");
            gr_mat_print(X, ctx);
            flint_printf("AX:\n");
            gr_mat_print(AX, ctx);
            flint_printf("\n");
            fflush(stdout);
            flint_abort();
        }

        if ((ctx->which_ring == GR_CTX_FMPQ || ctx->which_ring == GR_CTX_FQ ||
              ctx->which_ring == GR_CTX_FQ_NMOD || ctx->which_ring == GR_CTX_FQ_ZECH) && status == GR_UNABLE)
        {
            flint_printf("FAIL: unable over exact field\n");
            fflush(stdout);
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(X, ctx);
        gr_mat_clear(AX, ctx);

        gr_ctx_clear(ctx);
    }

    /* test random solvable systems */
    for (i = 0; i < 1000; i++)
    {
        status = GR_SUCCESS;

/*
        fmpz p = 17;
        gr_ctx_init_fq_nmod(ctx, &p, 1, "a");
*/

        while (1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_field(ctx) == T_TRUE)
                break;

            gr_ctx_clear(ctx);
        }

        m = n_randint(state, 6);
        n = n_randint(state, 6);
        k = n_randint(state, 6);

        gr_mat_init(A, m, k, ctx);
        gr_mat_init(B, m, n, ctx);
        gr_mat_init(X, k, n, ctx);
        gr_mat_init(X2, k, n, ctx);
        gr_mat_init(AX, m, n, ctx);

/*
        status |= gr_mat_randrank(A, state, n_randint(state, FLINT_MIN(m, k) + 1), ctx);
*/
        status |= gr_mat_randtest(A, state, ctx);

        status |= gr_mat_randtest(X2, state, ctx);
        status |= gr_mat_mul(B, A, X2, ctx);

        status |= gr_mat_solve_field(X, A, B, ctx);
        status |= gr_mat_mul(AX, A, X, ctx);

        if (status == GR_DOMAIN || (status == GR_SUCCESS && gr_mat_equal(AX, B, ctx) == T_FALSE))
        {
            flint_printf("FAIL:\n");
            flint_printf("status = %d\n", status);
            flint_printf("AX != B!\n");
            flint_printf("A:\n");
            gr_mat_print(A, ctx);
            flint_printf("B:\n");
            gr_mat_print(B, ctx);
            flint_printf("X:\n");
            gr_mat_print(X, ctx);
            flint_printf("AX:\n");
            gr_mat_print(AX, ctx);
            flint_printf("\n");

            fflush(stdout);
            flint_abort();
        }

        if ((ctx->which_ring == GR_CTX_FMPQ || ctx->which_ring == GR_CTX_FQ ||
              ctx->which_ring == GR_CTX_FQ_NMOD || ctx->which_ring == GR_CTX_FQ_ZECH) && status == GR_UNABLE)
        {
            flint_printf("FAIL: unable over exact field\n");
            fflush(stdout);
            flint_abort();
        }

        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_mat_clear(X, ctx);
        gr_mat_clear(X2, ctx);
        gr_mat_clear(AX, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
