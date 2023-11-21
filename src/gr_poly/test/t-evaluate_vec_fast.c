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
#include "gr_vec.h"
#include "gr_poly.h"

FLINT_DLL extern gr_static_method_table _ca_methods;

TEST_FUNCTION_START(gr_poly_evaluate_vec_fast, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_poly_t f;
        gr_vec_t x, y, z;
        slong n, m;

        gr_ctx_init_random(ctx, state);

        if (ctx->methods == _ca_methods)
        {
            n = n_randint(state, 8);
            m = n_randint(state, 8);
        }
        else
        {
            n = n_randint(state, 30);
            m = n_randint(state, 30);
        }

        gr_poly_init(f, ctx);
        gr_vec_init(x, n, ctx);
        gr_vec_init(y, n, ctx);
        gr_vec_init(z, n, ctx);

        status |= gr_poly_randtest(f, state, m, ctx);
        status |= _gr_vec_randtest(x->entries, state, n, ctx);
        status |= _gr_vec_randtest(y->entries, state, n, ctx);

        status |= gr_poly_evaluate_vec_fast(y, f, x, ctx);
        status |= gr_poly_evaluate_vec_iter(z, f, x, ctx);

        if (status == GR_SUCCESS && _gr_vec_equal(y->entries, z->entries, n, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n\n");
            flint_printf("f = "); gr_poly_print(f, ctx); flint_printf("\n");
            flint_printf("x = "); gr_vec_print(x, ctx); flint_printf("\n");
            flint_printf("y = "); gr_vec_print(y, ctx); flint_printf("\n");
            flint_printf("z = "); gr_vec_print(z, ctx); flint_printf("\n");
            flint_abort();
        }

        gr_poly_clear(f, ctx);
        gr_vec_clear(x, ctx);
        gr_vec_clear(y, ctx);
        gr_vec_clear(z, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
