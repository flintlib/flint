/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_gen, state)
{
    slong iter;

    for (iter = 0; iter < 100; iter++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f1, f2;
        slong len, k1, k2;
        flint_bitcnt_t exp_bits;
        int status;
        truth_t eq1, eq2, eq3, eq4;

        gr_ctx_init_random(cctx, state);
        gr_mpoly_ctx_init_rand(ctx, state, cctx, 10);

        gr_mpoly_init(f1, ctx);
        gr_mpoly_init(f2, ctx);

        len = n_randint(state, 10);
        exp_bits = n_randint(state, 100) + 2;

        status = gr_mpoly_zero(f1, ctx);
        status = gr_mpoly_one(f2, ctx);

        if (gr_mpoly_equal(f1, f2, ctx) != T_TRUE)
        {
            gr_mpoly_randtest_bits(f1, state, len, exp_bits, ctx);
            gr_mpoly_randtest_bits(f2, state, len, exp_bits, ctx);

            k1 = n_randint(state, FLINT_MAX(1, GR_MPOLY_NVARS(ctx)));
            k2 = n_randint(state, FLINT_MAX(1, GR_MPOLY_NVARS(ctx)));

            status = GR_SUCCESS;

            status |= gr_mpoly_gen(f1, k1, ctx);
            gr_mpoly_assert_canonical(f1, ctx);

            status |= gr_mpoly_gen(f2, k2, ctx);
            gr_mpoly_assert_canonical(f2, ctx);

            eq1 = gr_mpoly_is_gen(f1, k1, ctx);
            eq2 = gr_mpoly_is_gen(f1, -1, ctx);
            eq3 = gr_mpoly_is_gen(f2, k2, ctx);
            eq4 = gr_mpoly_is_gen(f2, -1, ctx);

            if (status == GR_SUCCESS && (eq1 == T_FALSE || eq2 == T_FALSE || eq3 == T_FALSE || eq4 == T_FALSE))
            {
                flint_printf("FAIL\n");
                gr_ctx_println(cctx);
                flint_printf("nvars = %wd\n", GR_MPOLY_NVARS(ctx));
                flint_printf("f1 = "); gr_mpoly_print_pretty(f1, ctx); flint_printf("\n");
                flint_printf("f2 = "); gr_mpoly_print_pretty(f2, ctx); flint_printf("\n");
                flint_printf("k1 = %wd\n", k1);
                flint_printf("k2 = %wd\n", k2);
                flint_printf("eq = %d %d %d %d\n", eq1, eq2, eq3, eq4);
                fflush(stdout);
                flint_abort();
            }

            if (n_randint(state, 2))
                status |= gr_mpoly_add(f1, f1, f2, ctx);
            else
                status |= gr_mpoly_mul(f1, f1, f2, ctx);

            eq1 = gr_mpoly_is_gen(f1, k1, ctx);
            eq2 = gr_mpoly_is_gen(f1, k2, ctx);
            eq3 = gr_mpoly_is_gen(f1, -1, ctx);

            if (status == GR_SUCCESS && (eq1 == T_TRUE || eq2 == T_TRUE || eq3 == T_TRUE))
            {
                flint_printf("FAIL\n");
                gr_ctx_println(cctx);
                flint_printf("nvars = %wd\n", GR_MPOLY_NVARS(ctx));
                flint_printf("f1 = "); gr_mpoly_print_pretty(f1, ctx); flint_printf("\n");
                flint_printf("f2 = "); gr_mpoly_print_pretty(f2, ctx); flint_printf("\n");
                flint_printf("k1 = %wd\n", k1);
                flint_printf("k2 = %wd\n", k2);
                flint_printf("eq = %d %d %d\n", eq1, eq2, eq3);
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f1, ctx);
        gr_mpoly_clear(f2, ctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
