/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "gr_vec.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_interpolate_newton, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_vec_t X, Y, Z;
        gr_poly_t F, G, H;
        slong n;

        gr_ctx_init_random(ctx, state);
        while (gr_ctx_is_integral_domain(ctx) != T_TRUE || ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        gr_poly_init(F, ctx);
        gr_poly_init(G, ctx);
        gr_poly_init(H, ctx);

        status |= gr_poly_randtest(F, state, n_randint(state, 6), ctx);
        status |= gr_poly_randtest(G, state, n_randint(state, 6), ctx);
        status |= gr_poly_randtest(H, state, n_randint(state, 6), ctx);

        n = F->length + n_randint(state, 3);

        gr_vec_init(X, n, ctx);
        gr_vec_init(Y, n, ctx);
        gr_vec_init(Z, n, ctx);

        /* Test recovering random polynomial */
        status |= _gr_vec_randtest(X->entries, state, n, ctx);
        status |= gr_poly_evaluate_vec_iter(Y, F, X, ctx);

        status |= gr_poly_interpolate_newton(G, X, Y, ctx);

        if (status == GR_SUCCESS)
        {
            if (gr_poly_equal(F, G, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                gr_ctx_println(ctx);
                flint_printf("X = "); gr_vec_print(X, ctx); flint_printf("\n");
                flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
                flint_printf("Y = "); gr_vec_print(Y, ctx); flint_printf("\n");
                flint_printf("G = "); gr_poly_print(G, ctx); flint_printf("\n");
                flint_abort();
            }

            status |= gr_poly_interpolate_exact_newton(H, X, Y, ctx);

            if (status != GR_SUCCESS || gr_poly_equal(F, H, ctx) == T_FALSE)
            {
                if (gr_poly_equal(F, H, ctx) == T_FALSE)
                {
                    flint_printf("FAIL (exact)\n\n");
                    gr_ctx_println(ctx);
                    flint_printf("X = "); gr_vec_print(X, ctx); flint_printf("\n");
                    flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
                    flint_printf("Y = "); gr_vec_print(Y, ctx); flint_printf("\n");
                    flint_printf("H = "); gr_poly_print(H, ctx); flint_printf("\n");
                    flint_abort();
                }
            }
        }

        /* Test interpolation from random points */
        status |= _gr_vec_randtest(X->entries, state, n, ctx);
        status |= _gr_vec_randtest(Y->entries, state, n, ctx);

        status |= gr_poly_interpolate_newton(G, X, Y, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_evaluate_vec_iter(Z, G, X, ctx);

            if (status == GR_SUCCESS && _gr_vec_equal(Z->entries, Y->entries, Y->length, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                gr_ctx_println(ctx);
                flint_printf("X = "); gr_vec_print(X, ctx); flint_printf("\n");
                flint_printf("Y = "); gr_vec_print(Y, ctx); flint_printf("\n");
                flint_printf("G = "); gr_poly_print(F, ctx); flint_printf("\n");
                flint_printf("Z = "); gr_vec_print(Z, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(F, ctx);
        gr_poly_clear(G, ctx);
        gr_poly_clear(H, ctx);
        gr_vec_clear(X, ctx);
        gr_vec_clear(Y, ctx);
        gr_vec_clear(Z, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
