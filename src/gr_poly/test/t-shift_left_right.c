/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_poly.h"

TEST_FUNCTION_START(gr_poly_shift_left_right, state)
{
    int i;

    /* Check aliasing of a and b for left shift */
    for (i = 0; i < 100; i++)
    {
        gr_ctx_t ctx;
        gr_poly_t a, b;
        int status = GR_SUCCESS;
        slong shift = n_randint(state, 5);

        gr_ctx_init_random(ctx, state);
        gr_poly_init(a, ctx);
        gr_poly_init(b, ctx);
        status |= gr_poly_randtest(a, state, n_randint(state, 5), ctx);
        status |= gr_poly_randtest(b, state, n_randint(state, 5), ctx);

        status |= gr_poly_shift_left(b, a, shift, ctx);
        status |= gr_poly_shift_left(a, a, shift, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(a, b, ctx) == T_FALSE)
        {
            flint_printf("FAIL (alias left):\n");
            gr_poly_print(a, ctx), flint_printf("\n\n");
            gr_poly_print(b, ctx), flint_printf("\n\n");
            flint_abort();
        }

        gr_poly_clear(a, ctx);
        gr_poly_clear(b, ctx);
        gr_ctx_clear(ctx);
    }

    /* Check aliasing of a and b for right shift */
    for (i = 0; i < 100; i++)
    {
        gr_ctx_t ctx;
        gr_poly_t a, b;
        int status = GR_SUCCESS;
        slong shift = n_randint(state, 5);

        gr_ctx_init_random(ctx, state);
        gr_poly_init(a, ctx);
        gr_poly_init(b, ctx);
        status |= gr_poly_randtest(a, state, n_randint(state, 5), ctx);
        status |= gr_poly_randtest(b, state, n_randint(state, 5), ctx);

        status |= gr_poly_shift_right(b, a, shift, ctx);
        status |= gr_poly_shift_right(a, a, shift, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(a, b, ctx) == T_FALSE)
        {
            flint_printf("FAIL (alias right):\n");
            gr_poly_print(a, ctx), flint_printf("\n\n");
            gr_poly_print(b, ctx), flint_printf("\n\n");
            flint_abort();
        }

        gr_poly_clear(a, ctx);
        gr_poly_clear(b, ctx);
        gr_ctx_clear(ctx);
    }

    /* Check shift left then right does nothing */
    for (i = 0; i < 100; i++)
    {
        gr_ctx_t ctx;
        gr_poly_t a, b, c;
        int status = GR_SUCCESS;
        slong shift = n_randint(state, 5);

        gr_ctx_init_random(ctx, state);
        gr_poly_init(a, ctx);
        gr_poly_init(b, ctx);
        gr_poly_init(c, ctx);
        status |= gr_poly_randtest(a, state, n_randint(state, 5), ctx);
        status |= gr_poly_randtest(b, state, n_randint(state, 5), ctx);
        status |= gr_poly_randtest(c, state, n_randint(state, 5), ctx);

        status |= gr_poly_shift_left(b, a, shift, ctx);
        status |= gr_poly_shift_right(c, b, shift, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(c, a, ctx) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            gr_poly_print(a, ctx), flint_printf("\n\n");
            gr_poly_print(b, ctx), flint_printf("\n\n");
            gr_poly_print(c, ctx), flint_printf("\n\n");
            flint_abort();
        }

        gr_poly_clear(a, ctx);
        gr_poly_clear(b, ctx);
        gr_poly_clear(c, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
