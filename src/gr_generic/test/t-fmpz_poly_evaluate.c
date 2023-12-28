/*
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_poly.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_generic.h"

static int
evaluate(flint_rand_t state, gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    switch (n_randint(state, 4))
    {
        case 0:
            status |= gr_fmpz_poly_evaluate(res, f, x, ctx);
            break;
        case 1:
            status |= gr_set(res, x, ctx);
            status |= gr_fmpz_poly_evaluate(res, f, res, ctx);
            break;
        case 2:
            status |= gr_fmpz_poly_evaluate_horner(res, f, x, ctx);
            break;
        case 3:
            status |= gr_set(res, x, ctx);
            status |= gr_fmpz_poly_evaluate_horner(res, f, res, ctx);
            break;
        case 4:
            status |= gr_fmpz_poly_evaluate_rectangular(res, f, x, ctx);
            break;
        default:
            status |= gr_set(res, x, ctx);
            status |= gr_fmpz_poly_evaluate_rectangular(res, f, res, ctx);
            break;
    }

    return status;
}

TEST_FUNCTION_START(gr_generic_fmpz_poly_evaluate, state)
{
    slong iter;
    slong count_success = 0, count_unable = 0, count_domain = 0;

    for (iter = 0; iter < 10000; iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        fmpz_poly_t f, g, h;
        gr_ptr x;
        gr_ptr fx, gx, hx, y;

        gr_ctx_init_random(ctx, state);

        fmpz_poly_init(f);
        fmpz_poly_init(g);
        fmpz_poly_init(h);

        x = gr_heap_init(ctx);
        fx = gr_heap_init(ctx);
        gx = gr_heap_init(ctx);
        hx = gr_heap_init(ctx);
        y = gr_heap_init(ctx);

        GR_MUST_SUCCEED(gr_randtest(x, state, ctx));

        fmpz_poly_randtest(f, state, 1 + n_randint(state, 30), 1 + n_randint(state, 100));
        fmpz_poly_randtest(g, state, 1 + n_randint(state, 30), 1 + n_randint(state, 100));
        fmpz_poly_add(h, f, g);

        status |= evaluate(state, fx, f, x, ctx);
        status |= evaluate(state, gx, g, x, ctx);
        status |= evaluate(state, hx, h, x, ctx);
        status |= gr_add(y, fx, gx, ctx);

        if (status == GR_SUCCESS && gr_equal(y, hx, ctx) == T_FALSE)
        {
            flint_printf("FAIL!\n");
            flint_printf("f = "); fmpz_poly_print(f); flint_printf("\n\n");
            flint_printf("g = "); fmpz_poly_print(g); flint_printf("\n\n");
            flint_printf("h = "); fmpz_poly_print(h); flint_printf("\n\n");
            flint_printf("x = "); gr_print(x, ctx); flint_printf("\n\n");
            flint_printf("fx = "); gr_print(fx, ctx); flint_printf("\n\n");
            flint_printf("gx = "); gr_print(gx, ctx); flint_printf("\n\n");
            flint_printf("hx = "); gr_print(hx, ctx); flint_printf("\n\n");
            flint_printf("y = "); gr_print(y, ctx); flint_printf("\n\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        fmpz_poly_clear(f);
        fmpz_poly_clear(g);
        fmpz_poly_clear(h);

        gr_heap_clear(x, ctx);
        gr_heap_clear(fx, ctx);
        gr_heap_clear(gx, ctx);
        gr_heap_clear(hx, ctx);
        gr_heap_clear(y, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_unable, count_domain);
}
