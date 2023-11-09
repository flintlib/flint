/*
    Copyright (C) 2023 Joel Dahne

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "flint.h"
#include "gr_poly.h"
#include "ulong_extras.h"

TEST_FUNCTION_START(gr_poly_nth_derivative, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        ulong nth;
        slong j;
        gr_ctx_t ctx;
        gr_poly_t a, b;

        gr_ctx_init_random(ctx, state);

        gr_poly_init(a, ctx);
        gr_poly_init(b, ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(a, state, n_randint(state, 30), ctx);
        status |= gr_poly_randtest(b, state, n_randint(state, 30), ctx);

        nth = n_randint(state, 30);

        if (n_randint(state, 2)) {
            status |= gr_poly_nth_derivative(b, a, nth, ctx);
        }
        else
        {
            status |= gr_poly_set(b, a, ctx);
            status |= gr_poly_nth_derivative(b, b, nth, ctx);
        }

        /* Compute derivative iteratively */
        for (j = 0; j < nth; j ++)
            status |= gr_poly_derivative(a, a, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(a, b, ctx) == T_FALSE)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n");
            flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(a, ctx);
        gr_poly_clear(b, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
