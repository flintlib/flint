/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("make_monic....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t f, fc, g, h;
        gr_ptr c;

        gr_ctx_init_random(ctx, state);

        gr_poly_init(f, ctx);
        gr_poly_init(fc, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(h, ctx);
        c = gr_heap_init(ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(f, state, 1 + n_randint(state, 6), ctx);
        status |= gr_poly_randtest(fc, state, 1 + n_randint(state, 6), ctx);
        status |= gr_poly_randtest(g, state, 1 + n_randint(state, 6), ctx);
        status |= gr_randtest(c, state, ctx);

        status |= gr_poly_mul_scalar(fc, f, c, ctx);

        if (n_randint(state, 2))
        {
            status |= gr_poly_make_monic(g, f, ctx);
        }
        else
        {
            status |= gr_poly_set(g, f, ctx);
            status |= gr_poly_make_monic(g, g, ctx);
        }

        if (status == GR_SUCCESS && gr_poly_is_monic(g, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
            flint_abort();
        }

        status |= gr_poly_make_monic(h, fc, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(g, h, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
            flint_printf("h = "); gr_poly_print(h, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(fc, ctx);
        gr_poly_clear(g, ctx);
        gr_poly_clear(h, ctx);
        gr_heap_clear(c, ctx);

        gr_ctx_clear(ctx);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
