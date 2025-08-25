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

TEST_FUNCTION_START(gr_poly_product_roots, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_vec_t R, Y;
        gr_poly_t F;
        slong n;

        gr_ctx_init_random(ctx, state);
        /* Hack: avoid because slow */
        while (ctx->methods == _ca_methods)
        {
            gr_ctx_clear(ctx);
            gr_ctx_init_random(ctx, state);
        }

        n = n_randint(state, 10);

        gr_vec_init(R, n, ctx);
        gr_vec_init(Y, n, ctx);
        gr_poly_init(F, ctx);

        status |= _gr_vec_randtest(R->entries, state, n, ctx);
        status |= gr_poly_randtest(F, state, n_randint(state, 6), ctx);
        status |= gr_poly_product_roots(F, R, ctx);
        status |= gr_poly_evaluate_vec_iter(Y, F, R, ctx);

        if (status == GR_SUCCESS && (_gr_vec_is_zero(Y->entries, n, ctx) == T_FALSE ||
            (gr_ctx_is_zero_ring(ctx) == T_FALSE && F->length != n + 1)))
        {
            flint_printf("FAIL\n\n");
            gr_ctx_println(ctx);
            flint_printf("R = "); gr_vec_print(R, ctx); flint_printf("\n");
            flint_printf("F = "); gr_poly_print(F, ctx); flint_printf("\n");
            flint_printf("Y = "); gr_vec_print(Y, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_vec_clear(R, ctx);
        gr_vec_clear(Y, ctx);
        gr_poly_clear(F, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
