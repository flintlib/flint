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

FLINT_DLL extern gr_static_method_table _ca_methods;

int
test_inv_series(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    slong n;
    gr_poly_t A, B, AB, one;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(AB, ctx);
    gr_poly_init(one, ctx);

    if (ctx->methods == _ca_methods)
        n = n_randint(state, 5);
    else
        n = n_randint(state, 20);

    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 20, ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(B, state, 20, ctx));

    switch (which)
    {
        case 0:
            status |= gr_poly_inv_series(B, A, n, ctx);
            break;
        case 1:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_inv_series(B, B, n, ctx);
            break;
        case 2:
            status |= gr_poly_inv_series_basecase(B, A, n, ctx);
            break;
        case 3:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_inv_series_basecase(B, B, n, ctx);
            break;
        case 4:
            status |= gr_poly_inv_series_newton(B, A, n, n_randint(state, 20), ctx);
            break;
        case 5:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_inv_series_newton(B, B, n, n_randint(state, 20), ctx);
            break;
        default:
            flint_abort();
    }

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_mullow(AB, A, B, n, ctx);
        status |= gr_poly_one(one, ctx);
        status |= gr_poly_truncate(one, one, n, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(AB, one, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("which = %d\n\n", which);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("AB = "); gr_poly_print(AB, ctx); flint_printf("\n");
            flint_abort();
        }
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(B, ctx);
    gr_poly_clear(AB, ctx);
    gr_poly_clear(one, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_inv_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test_inv_series(state, 0);
        test_inv_series(state, 1);
        test_inv_series(state, 2);
        test_inv_series(state, 3);
        test_inv_series(state, 4);
        test_inv_series(state, 5);
    }

    TEST_FUNCTION_END(state);
}
