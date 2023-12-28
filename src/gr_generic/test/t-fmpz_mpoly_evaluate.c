/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_generic.h"

TEST_FUNCTION_START(gr_generic_fmpz_mpoly_evaluate, state)
{
    slong iter;
    slong count_success = 0, count_unable = 0, count_domain = 0;

    for (iter = 0; iter < 10000; iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        slong i, n;
        fmpz_mpoly_ctx_t mctx;
        fmpz_mpoly_t f, g, h;
        gr_ptr x;
        gr_ptr fx, gx, hx, y;
        slong sz;

        gr_ctx_init_random(ctx, state);
        sz = ctx->sizeof_elem;

        n = 1 + n_randint(state, 5);
        fmpz_mpoly_ctx_init(mctx, n, ORD_LEX);

        fmpz_mpoly_init(f, mctx);
        fmpz_mpoly_init(g, mctx);
        fmpz_mpoly_init(h, mctx);

        x = gr_heap_init_vec(n, ctx);
        fx = gr_heap_init(ctx);
        gx = gr_heap_init(ctx);
        hx = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        status |= _gr_vec_randtest(x, state, n, ctx);

        fmpz_mpoly_randtest_bound(f, state, 1 + n_randint(state, 30), 10, 1 + n_randint(state, 6), mctx);
        fmpz_mpoly_randtest_bound(g, state, 1 + n_randint(state, 30), 10, 1 + n_randint(state, 6), mctx);
        fmpz_mpoly_add(h, f, g, mctx);

        status |= gr_fmpz_mpoly_evaluate(fx, f, x, mctx, ctx);
        status |= gr_fmpz_mpoly_evaluate(gx, g, x, mctx, ctx);
        status |= gr_fmpz_mpoly_evaluate(hx, h, x, mctx, ctx);
        status |= gr_add(y, fx, gx, ctx);

        if (status == GR_SUCCESS && gr_equal(y, hx, ctx) == T_FALSE)
        {
            flint_printf("FAIL!\n");
            flint_printf("f = "); fmpz_mpoly_print_pretty(f, NULL, mctx); flint_printf("\n\n");
            flint_printf("g = "); fmpz_mpoly_print_pretty(g, NULL, mctx); flint_printf("\n\n");
            flint_printf("h = "); fmpz_mpoly_print_pretty(h, NULL, mctx); flint_printf("\n\n");
            for (i = 0; i < n; i++)
            {
                flint_printf("x%wd = ", i + 1); gr_print(GR_ENTRY(x, i, sz), ctx); flint_printf("\n\n");
            }
            flint_printf("fx = "); gr_print(fx, ctx); flint_printf("\n\n");
            flint_printf("gx = "); gr_print(gx, ctx); flint_printf("\n\n");
            flint_printf("hx = "); gr_print(hx, ctx); flint_printf("\n\n");
            flint_printf("y = "); gr_print(y, ctx); flint_printf("\n\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        fmpz_mpoly_clear(f, mctx);
        fmpz_mpoly_clear(g, mctx);
        fmpz_mpoly_clear(h, mctx);
        fmpz_mpoly_ctx_clear(mctx);

        gr_heap_clear_vec(x, n, ctx);
        gr_heap_clear(fx, ctx);
        gr_heap_clear(gx, ctx);
        gr_heap_clear(hx, ctx);
        gr_heap_clear(y, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
