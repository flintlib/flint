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

TEST_FUNCTION_START(gr_poly_div_root, state)
{
    slong iter;

    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        int status;
        gr_ctx_t ctx;
        gr_poly_t A, B, Q, Q2, R2;
        gr_ptr c, R;

        gr_ctx_init_random_commutative_ring(ctx, state);

        gr_poly_init(A, ctx);
        gr_poly_init(B, ctx);
        gr_poly_init(Q, ctx);
        gr_poly_init(Q2, ctx);
        gr_poly_init(R2, ctx);
        c = gr_heap_init(ctx);
        R = gr_heap_init(ctx);

        status = GR_SUCCESS;

        status |= gr_poly_randtest(A, state, 1 + n_randint(state, 6), ctx);
        status |= gr_randtest(c, state, ctx);
        status |= gr_randtest(R, state, ctx);

        /* B = x - c */
        status |= gr_poly_set_coeff_si(B, 1, 1, ctx);
        status |= gr_neg(c, c, ctx);
        status |= gr_poly_set_coeff_scalar(B, 0, c, ctx);
        status |= gr_neg(c, c, ctx);

        if (n_randint(state, 2))
        {
            status |= gr_poly_set(Q, A, ctx);
            status |= gr_poly_div_root(Q, R, Q, c, ctx);
        }
        else
        {
            status |= gr_poly_div_root(Q, R, A, c, ctx);
        }

        status |= gr_poly_divrem(Q2, R2, A, B, ctx);

        if (status == GR_SUCCESS)
        {
            status |= gr_poly_set_scalar(B, R, ctx);

            if (gr_poly_equal(Q, Q2, ctx) == T_FALSE || gr_poly_equal(B, R2, ctx) == T_FALSE)
            {
                flint_printf("FAIL\n\n");
                flint_printf("A = "); gr_poly_print(A, ctx); flint_printf("\n");
                flint_printf("c = "); gr_print(c, ctx); flint_printf("\n");
                flint_printf("Q = "); gr_poly_print(Q, ctx); flint_printf("\n");
                flint_printf("R = "); gr_print(R, ctx); flint_printf("\n");
                flint_printf("Q2 = "); gr_poly_print(Q2, ctx); flint_printf("\n");
                flint_printf("R2 = "); gr_poly_print(R2, ctx); flint_printf("\n");
                flint_abort();
            }
        }

        gr_poly_clear(A, ctx);
        gr_poly_clear(B, ctx);
        gr_poly_clear(Q, ctx);
        gr_poly_clear(Q2, ctx);
        gr_poly_clear(R2, ctx);
        gr_heap_clear(c, ctx);
        gr_heap_clear(R, ctx);

        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
