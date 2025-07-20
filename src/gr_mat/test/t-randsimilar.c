/*
    Copyright (C) 2025 Ricardo Buring

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_poly.h"
#include "gr_mat.h"

TEST_FUNCTION_START(gr_mat_randsimilar, state)
{
    slong iter;

    /* Check that random similarity transformations preserve the characteristic polynomial. */
    for (iter = 0; iter < 10000; iter++)
    {
        gr_ctx_t ctx;
        gr_mat_t A;
        gr_poly_t f, g;
        slong n, rank, randops_count, randsimilar_opcount;
        int status = GR_SUCCESS;

        gr_ctx_init_random(ctx, state);

        if (gr_ctx_is_integral_domain(ctx) == T_TRUE)
        {
            n = n_randint(state, 6);
            rank = n_randint(state, n + 1);
            randops_count = n_randint(state, 8);
            randsimilar_opcount = n_randint(state, 8);

            gr_mat_init(A, n, n, ctx);
            gr_poly_init(f, ctx);
            gr_poly_init(g, ctx);

            status |= gr_mat_randrank(A, state, rank, ctx);
            status |= gr_mat_randops(A, state, randops_count, ctx);
            status |= gr_mat_charpoly(f, A, ctx);
            status |= gr_mat_randsimilar(A, state, randsimilar_opcount, ctx);
            status |= gr_mat_charpoly(g, A, ctx);

            if (gr_poly_equal(f, g, ctx) == T_FALSE)
            {
                flint_printf("FAIL:\n");
                gr_ctx_println(ctx);
                flint_printf("A: "); gr_mat_print(A, ctx); flint_printf("\n");
                flint_printf("f: "); gr_poly_print(f, ctx); flint_printf("\n");
                flint_printf("g: "); gr_poly_print(g, ctx); flint_printf("\n");
                flint_abort();
            }

            gr_poly_clear(g, ctx);
            gr_poly_clear(f, ctx);
            gr_mat_clear(A, ctx);
        }

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
