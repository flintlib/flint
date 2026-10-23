/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_deflate_inflate, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t f, g, h, r;
        ulong n;
        slong i;

        gr_ctx_init_random(ctx, state);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);
        gr_poly_init(h, ctx);
        gr_poly_init(r, ctx);

        status = GR_SUCCESS;

        /* lengths on both sides of the algorithm cutoff in gr_poly_inflate */
        status |= gr_poly_randtest(f, state, n_randint(state, 2) ? n_randint(state, 8) : n_randint(state, 60), ctx);
        n = 1 + n_randint(state, 5);

        /* output with existing (possibly heap-allocated) entries */
        status |= gr_poly_randtest(g, state, n_randint(state, 200), ctx);
        status |= gr_poly_randtest(h, state, n_randint(state, 200), ctx);

        /* reference: f(x^n) coefficient by coefficient */
        for (i = 0; i < f->length; i++)
            status |= gr_poly_set_coeff_scalar(r, i * n, gr_poly_coeff_srcptr(f, i, ctx), ctx);

        if (n_randint(state, 2))
        {
            status |= gr_poly_inflate(g, f, n, ctx);
        }
        else
        {
            status |= gr_poly_set(g, f, ctx);
            status |= gr_poly_inflate(g, g, n, ctx);
        }

        if (status == GR_SUCCESS && gr_poly_equal(g, r, ctx) == T_FALSE)
        {
            flint_printf("FAIL (inflate)\n\n");
            gr_ctx_println(ctx); flint_printf("n = %wu\n", n);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
            flint_printf("r = "); gr_poly_print(r, ctx); flint_printf("\n");
            flint_abort();
        }

        if (n_randint(state, 2))
        {
            status |= gr_poly_deflate(h, g, n, ctx);
        }
        else
        {
            status |= gr_poly_set(h, g, ctx);
            status |= gr_poly_deflate(h, h, n, ctx);
        }

        if (status == GR_SUCCESS && gr_poly_equal(h, f, ctx) == T_FALSE)
        {
            flint_printf("FAIL (deflate)\n\n");
            gr_ctx_println(ctx); flint_printf("n = %wu\n", n);
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
            flint_printf("g = "); gr_poly_print(g, ctx); flint_printf("\n");
            flint_printf("h = "); gr_poly_print(h, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_poly_clear(h, ctx);
        gr_poly_clear(r, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
