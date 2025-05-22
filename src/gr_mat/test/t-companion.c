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

TEST_FUNCTION_START(gr_mat_companion, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        int status = GR_SUCCESS;
        gr_ctx_t ctx;
        gr_mat_t A, B;
        gr_poly_t f, g;
        gr_ptr den;
        slong n;

        do {
            gr_ctx_init_random(ctx, state);
        } while (gr_ctx_is_field(ctx) != T_TRUE);

        gr_poly_init(f, ctx);
        gr_poly_init(g, ctx);

        do {
            status |= gr_poly_randtest(f, state, 1 + n_randint(state, 8), ctx);
            n = gr_poly_length(f, ctx) - 1;
        } while (n < 0 || gr_is_invertible(gr_poly_coeff_srcptr(f, n, ctx), ctx) != T_TRUE);

        gr_mat_init(A, n, n, ctx);
        status |= gr_mat_randtest(A, state, ctx);
        status |= gr_mat_companion(A, f, ctx);
        status |= gr_mat_charpoly(g, A, ctx);

        status |= gr_poly_scalar_mul(g, gr_poly_coeff_srcptr(f, n, ctx), g, ctx);

        if (gr_poly_equal(g, f, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            flint_printf("A:\n"), gr_mat_print(A, ctx), flint_printf("\n");
            flint_printf("f:\n"), gr_poly_print(f, ctx), flint_printf("\n");
            flint_printf("g:\n"), gr_poly_print(g, ctx), flint_printf("\n");
            flint_abort();
        }

        GR_TMP_INIT(den, ctx);
        gr_mat_init(B, n, n, ctx);
        status |= gr_mat_randtest(B, state, ctx);
        status |= gr_mat_companion_fraction(B, den, f, ctx);
        status |= gr_mat_mul_scalar(A, A, den, ctx);

        if (gr_mat_equal(A, B, ctx) == T_FALSE)
        {
            flint_printf("FAIL\n");
            flint_printf("A:\n"), gr_mat_print(A, ctx), flint_printf("\n");
            flint_printf("B:\n"), gr_mat_print(B, ctx), flint_printf("\n");
            flint_printf("f:\n"), gr_poly_print(f, ctx), flint_printf("\n");
            flint_printf("g:\n"), gr_poly_print(g, ctx), flint_printf("\n");
            flint_abort();
        }

        GR_TMP_CLEAR(den, ctx);
        gr_mat_clear(A, ctx);
        gr_mat_clear(B, ctx);
        gr_poly_clear(f, ctx);
        gr_poly_clear(g, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
