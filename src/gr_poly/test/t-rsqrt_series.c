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
test_rsqrt_series(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    slong n;
    gr_poly_t A, B, C;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(C, ctx);

    if (ctx->methods == _ca_methods)
        n = n_randint(state, 5);
    else
        n = n_randint(state, 20);

    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 20, ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(B, state, 20, ctx));

    if (n_randint(state, 2))
    {
        status |= gr_poly_mullow(A, A, A, n_randint(state, 20), ctx);
        status |= gr_poly_inv_series(A, A, n_randint(state, 20), ctx);
    }

    switch (which)
    {
        case 0:
            status |= gr_poly_rsqrt_series(B, A, n, ctx);
            break;
        case 1:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_rsqrt_series(B, B, n, ctx);
            break;
        case 2:
            status |= gr_poly_rsqrt_series_basecase(B, A, n, ctx);
            break;
        case 3:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_rsqrt_series_basecase(B, B, n, ctx);
            break;
        case 4:
            status |= gr_poly_rsqrt_series_newton(B, A, n, n_randint(state, 20), ctx);
            break;
        case 5:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_rsqrt_series_newton(B, B, n, n_randint(state, 20), ctx);
            break;
        case 6:
            status |= gr_poly_rsqrt_series_miller(B, A, n, ctx);
            break;
        case 7:
            status |= gr_poly_set(B, A, ctx);
            status |= gr_poly_rsqrt_series_miller(B, B, n, ctx);
            break;
        default:
            flint_abort();
    }

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_mullow(C, B, B, n, ctx);
        status |= gr_poly_inv_series(C, C, n, ctx);
        status |= gr_poly_truncate(A, A, n, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(C, A, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("which = %d\n\n", which);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
            flint_abort();
        }
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(B, ctx);
    gr_poly_clear(C, ctx);

    gr_ctx_clear(ctx);

    return status;
}

TEST_FUNCTION_START(gr_poly_rsqrt_series, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        test_rsqrt_series(state, n_randint(state, 8));
    }

    TEST_FUNCTION_END(state);
}
