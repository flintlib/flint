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

TEST_FUNCTION_START(gr_poly_mulmid_classical, state)
{
    slong iter;

    for (iter = 0; iter < 1000; iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t a, b, res1, res2;
        slong nlo, nhi, N;
        int aliasing;

        gr_ctx_init_random(ctx, state);
        if (gr_ctx_is_finite(ctx) == T_TRUE || gr_ctx_has_real_prec(ctx) == T_TRUE)
            N = 10;
        else
            N = 4;

        gr_poly_init(a, ctx);
        gr_poly_init(b, ctx);
        gr_poly_init(res1, ctx);
        gr_poly_init(res2, ctx);

        status = GR_SUCCESS;

        nlo = n_randint(state, N);
        nhi = n_randint(state, N);

        status |= gr_poly_randtest(a, state, 1 + n_randint(state, N), ctx);
        status |= gr_poly_randtest(b, state, 1 + n_randint(state, N), ctx);
        status |= gr_poly_randtest(res1, state, 1 + n_randint(state, N), ctx);
        status |= gr_poly_randtest(res2, state, 1 + n_randint(state, N), ctx);
        aliasing = n_randint(state, 5);

        if (aliasing == 0)
        {
            status |= gr_poly_mulmid_classical(res1, a, b, nlo, nhi, ctx);
        }
        else if (aliasing == 1)
        {
            status |= gr_poly_set(res1, a, ctx);
            status |= gr_poly_mulmid_classical(res1, res1, b, nlo, nhi, ctx);
        }
        else if (aliasing == 2)
        {
            status |= gr_poly_set(res1, b, ctx);
            status |= gr_poly_mulmid_classical(res1, a, res1, nlo, nhi, ctx);
        }
        else if (aliasing == 3)
        {
            status |= gr_poly_mulmid_classical(res1, a, a, nlo, nhi, ctx);
        }
        else if (aliasing == 4)
        {
            status |= gr_poly_set(res1, a, ctx);
            status |= gr_poly_mulmid_classical(res1, res1, res1, nlo, nhi, ctx);
        }

        if (status == GR_SUCCESS)
        {
            if (aliasing == 3 || aliasing == 4)
                status |= gr_poly_mullow(res2, a, a, nhi, ctx);
            else
                status |= gr_poly_mullow(res2, a, b, nhi, ctx);

            status |= gr_poly_shift_right(res2, res2, nlo, ctx);

            if (status == GR_SUCCESS && gr_poly_equal(res1, res2, ctx) == T_FALSE)
            {
                flint_printf("FAIL: gr_poly_mulmid\n\n");
                gr_ctx_println(ctx);
                flint_printf("nlo = %wd, nhi = %wd, aliasing = %d\n", nlo, nhi, aliasing);
                flint_printf("a = "); gr_poly_print(a, ctx); flint_printf("\n");
                flint_printf("b = "); gr_poly_print(b, ctx); flint_printf("\n");
                flint_printf("res1 = "); gr_poly_print(res1, ctx); flint_printf("\n");
                flint_printf("res2 = "); gr_poly_print(res2, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(a, ctx);
        gr_poly_clear(b, ctx);
        gr_poly_clear(res1, ctx);
        gr_poly_clear(res2, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
