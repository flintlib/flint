/*
    Copyright (C) 2010, 2022 Fredrik Johansson
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2015 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mat.h"
#include "gr_poly.h"
/* #include "fq.h" */

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_GR_FUNCTION_START(gr_mat_minpoly_field, state, count_success, count_unable, count_domain)
{
    slong iter;

    /* minpoly(A) divides charpoly(A) */
    for (iter = 0; iter < 1000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        gr_poly_t p1, p2, q, r;
        slong m, n;
        int status = GR_SUCCESS;

/*
        while (1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_field(ctx) == T_TRUE)
                break;
            gr_ctx_clear(ctx);
        }
*/
        gr_ctx_init_random(ctx, state);

        if (ctx->methods == _ca_methods)
        {
            m = n_randint(state, 4);
            n = m;
        }
        else
        {
            m = n_randint(state, 10);
            n = m;
        }

        gr_poly_init(p1, ctx);
        gr_poly_init(p2, ctx);
        gr_poly_init(q, ctx);
        gr_poly_init(r, ctx);
        gr_mat_init(A, m, n, ctx);
        status |= gr_mat_randtest(A, state, ctx);

        status |= gr_mat_charpoly(p1, A, ctx);
        status |= gr_mat_minpoly_field(p2, A, ctx);

        if ((gr_ctx_is_field(ctx) == T_TRUE && status == GR_DOMAIN) || (ctx->which_ring == GR_CTX_FMPQ && status != GR_SUCCESS))
        {
            flint_printf("FAIL:\n");
            flint_printf("expected success\n");
            gr_ctx_println(ctx);
            gr_mat_print(A, ctx);
            fflush(stdout);
            flint_abort();
        }

        status |= gr_poly_divrem(q, r, p1, p2, ctx);

        if (status == GR_SUCCESS && gr_poly_is_zero(r, ctx) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            flint_printf("minpoly(A) doesn't divide charpoly(A).\n");
            fflush(stdout);
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        gr_mat_clear(A, ctx);
        gr_poly_clear(p1, ctx);
        gr_poly_clear(p2, ctx);
        gr_poly_clear(q, ctx);
        gr_poly_clear(r, ctx);

        gr_ctx_clear(ctx);
    }

    /* minpoly(P^{-1}AP) == minpoly(A) */
    for (iter = 0; iter < 1000; iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_mat_t A, B;
        gr_poly_t p1, p2;
        gr_ptr t;
        slong j, k, m, n;

        gr_ctx_init_random(ctx, state);

        if (ctx->methods == _ca_methods)
        {
            m = n_randint(state, 4);
            n = m;
        }
        else
        {
            m = n_randint(state, 10);
            n = m;
        }

        t = gr_heap_init(ctx);
        gr_poly_init(p1, ctx);
        gr_poly_init(p2, ctx);
        gr_mat_init(A, m, n, ctx);
        gr_mat_init(B, m, n, ctx);

        status |= gr_mat_randtest(A, state, ctx);

        for (j = 0; j < n / 2; j++)
        {
            for (k = 0; k < n / 2; k++)
            {
                status |= gr_zero(gr_mat_entry_ptr(A, j + n / 2, k, ctx), ctx);
                status |= gr_zero(gr_mat_entry_ptr(A, j, k + n / 2, ctx), ctx);
                status |= gr_set(gr_mat_entry_ptr(A, j + n / 2, k + n / 2, ctx), gr_mat_entry_ptr(A, j, k, ctx), ctx);
            }
        }

        status |= gr_mat_set(B, A, ctx);
        status |= gr_mat_minpoly_field(p1, A, ctx);

        if (status == GR_SUCCESS)
        {
            for (j = 0; j < n; j++)
            {
                status |= gr_set_si(t, n_randint(state, 6) - 3, ctx);
                status |= gr_mat_apply_row_similarity(B, n_randint(state, n), t, ctx);
            }

            status |= gr_mat_minpoly_field(p2, B, ctx);

            if (status == GR_SUCCESS && gr_poly_equal(p1, p2, ctx) == T_FALSE)
            {
                flint_printf("FAIL:\n");
                flint_printf("minpoly(P^{-1}AP) != minpoly(A).\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_heap_clear(t, ctx);
        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_poly_clear(p1, ctx);
        gr_poly_clear(p2, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
