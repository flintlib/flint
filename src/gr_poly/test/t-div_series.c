/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "gr_poly.h"

extern gr_static_method_table _ca_methods;

int
test_div_series(flint_rand_t state, int which)
{
    gr_ctx_t ctx;
    slong n;
    gr_poly_t A, B, C, D, E;
    int status = GR_SUCCESS;

    gr_ctx_init_random(ctx, state);

    gr_poly_init(A, ctx);
    gr_poly_init(B, ctx);
    gr_poly_init(C, ctx);
    gr_poly_init(D, ctx);
    gr_poly_init(E, ctx);

    if (ctx->methods == _ca_methods)
        n = n_randint(state, 5);
    else
        n = n_randint(state, 20);

    GR_MUST_SUCCEED(gr_poly_randtest(A, state, 1 + n_randint(state, 20), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(B, state, 1 + n_randint(state, 20), ctx));
    GR_MUST_SUCCEED(gr_poly_randtest(C, state, 1 + n_randint(state, 20), ctx));

    /* todo: randomly make exact multiple? */

    switch (which)
    {
        case 0:
            status |= gr_poly_div_series_basecase(C, A, B, n, ctx);
            break;
        case 1:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series_basecase(C, C, B, n, ctx);
            break;
        case 2:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series_basecase(C, A, C, n, ctx);
            break;

        case 3:
            status |= gr_poly_div_series_newton(C, A, B, n, n_randint(state, 20), ctx);
            break;
        case 4:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series_newton(C, C, B, n, n_randint(state, 20), ctx);
            break;
        case 5:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series_newton(C, A, C, n, n_randint(state, 20), ctx);
            break;

        case 6:
            status |= gr_poly_div_series(C, A, B, n, ctx);
            break;
        case 7:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series(C, C, B, n, ctx);
            break;
        case 8:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series(C, A, C, n, ctx);
            break;

        case 9:
            status |= gr_poly_div_series_invmul(C, A, B, n, ctx);
            break;
        case 10:
            status |= gr_poly_set(C, A, ctx);
            status |= gr_poly_div_series_invmul(C, C, B, n, ctx);
            break;
        case 11:
            status |= gr_poly_set(C, B, ctx);
            status |= gr_poly_div_series_invmul(C, A, C, n, ctx);
            break;


        default:
            abort();
    }

    if (status == GR_SUCCESS)
    {
        status |= gr_poly_mullow(D, C, B, n, ctx);
        status |= gr_poly_set(E, A, ctx);
        status |= gr_poly_truncate(E, n, ctx);

        if (status == GR_SUCCESS && gr_poly_equal(D, E, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("which = %d\n\n", which);
            gr_ctx_println(ctx);
            flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
            flint_printf("B = "); gr_poly_print(B, ctx); flint_printf("\n");
            flint_printf("C = "); gr_poly_print(C, ctx); flint_printf("\n");
            flint_printf("D = "); gr_poly_print(D, ctx); flint_printf("\n");
            flint_printf("E = "); gr_poly_print(E, ctx); flint_printf("\n");
            flint_abort();
        }
    }

    gr_poly_clear(A, ctx);
    gr_poly_clear(B, ctx);
    gr_poly_clear(C, ctx);
    gr_poly_clear(D, ctx);
    gr_poly_clear(E, ctx);

    gr_ctx_clear(ctx);

    return status;
}

int main()
{
    slong iter;
    flint_rand_t state;

    flint_printf("div_series....");
    fflush(stdout);

    flint_randinit(state);

    for (iter = 0; iter < 10000; iter++)
    {
        test_div_series(state, n_randint(state, 12));
    }

    flint_randclear(state);
    flint_cleanup_master();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
