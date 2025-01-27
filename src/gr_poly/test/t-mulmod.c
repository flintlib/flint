/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2011 Sebastian Pancratz
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_mulmod, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t a, b, res1, res2, f;

        while (1)
        {
            gr_ctx_init_random(ctx, state);
            if (gr_ctx_is_finite(ctx) == T_TRUE || gr_ctx_has_real_prec(ctx) == T_TRUE)
                break;
            gr_ctx_clear(ctx);
        }

        gr_poly_init(a, ctx);
        gr_poly_init(b, ctx);
        gr_poly_init(res1, ctx);
        gr_poly_init(res2, ctx);
        gr_poly_init(f, ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(a, state, 1 + n_randint(state, 20), ctx);
        status |= gr_poly_randtest(b, state, 1 + n_randint(state, 20), ctx);
        status |= gr_poly_randtest(res1, state, 1 + n_randint(state, 20), ctx);
        status |= gr_poly_randtest(res2, state, 1 + n_randint(state, 20), ctx);
        status |= gr_poly_randtest(f, state, 1 + n_randint(state, 20), ctx);

        /* test aliasing */
        switch (n_randint(state, 4))
        {
            case 0:
                status |= gr_poly_set(res1, a, ctx);
                status |= gr_poly_mulmod(res1, res1, b, f, ctx);
                break;
            case 1:
                status |= gr_poly_set(res1, b, ctx);
                status |= gr_poly_mulmod(res1, a, res1, f, ctx);
                break;
            case 2:
                status |= gr_poly_set(res1, f, ctx);
                status |= gr_poly_mulmod(res1, a, b, res1, ctx);
                break;
            default:
                status |= gr_poly_mulmod(res1, a, b, f, ctx);
                break;
        }

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_mul(res2, a, b, ctx);
            status |= gr_poly_rem(res2, res2, f, ctx);

            if (status == GR_SUCCESS && gr_poly_equal(res1, res2, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n");
                flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n");
                flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
                flint_printf("res1 = "); gr_poly_print(res1, ctx); flint_printf("\n");
                flint_printf("res2 = "); gr_poly_print(res2, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(a, ctx);
        gr_poly_clear(b, ctx);
        gr_poly_clear(res1, ctx);
        gr_poly_clear(res2, ctx);
        gr_poly_clear(f, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
