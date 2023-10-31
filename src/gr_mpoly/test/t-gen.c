/*
    Copyright (C) 2020 Daniel Schultz
    Copyright (C) 2022 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
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
        mpoly_ctx_t mctx;
        gr_mpoly_t f1, f2;
        slong len, k1, k2;
        flint_bitcnt_t exp_bits;
        int status;
        truth_t eq1, eq2, eq3, eq4;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 10);

        gr_mpoly_init(f1, mctx, cctx);
        gr_mpoly_init(f2, mctx, cctx);

        len = n_randint(state, 10);
        exp_bits = n_randint(state, 100) + 2;

        status = gr_mpoly_zero(f1, mctx, cctx);
        status = gr_mpoly_one(f2, mctx, cctx);

        if (gr_mpoly_equal(f1, f2, mctx, cctx) != T_TRUE)
        {
            gr_mpoly_randtest_bits(f1, state, len, exp_bits, mctx, cctx);
            gr_mpoly_randtest_bits(f2, state, len, exp_bits, mctx, cctx);

            k1 = n_randint(state, FLINT_MAX(1, mctx->nvars));
            k2 = n_randint(state, FLINT_MAX(1, mctx->nvars));

            status = GR_SUCCESS;

            status |= gr_mpoly_gen(f1, k1, mctx, cctx);
            gr_mpoly_assert_canonical(f1, mctx, cctx);

            status |= gr_mpoly_gen(f2, k2, mctx, cctx);
            gr_mpoly_assert_canonical(f2, mctx, cctx);

            eq1 = gr_mpoly_is_gen(f1, k1, mctx, cctx);
            eq2 = gr_mpoly_is_gen(f1, -1, mctx, cctx);
            eq3 = gr_mpoly_is_gen(f2, k2, mctx, cctx);
            eq4 = gr_mpoly_is_gen(f2, -1, mctx, cctx);

            if (status == GR_SUCCESS && (eq1 == T_FALSE || eq2 == T_FALSE || eq3 == T_FALSE || eq4 == T_FALSE))
            {
                flint_printf("FAIL\n");
                gr_ctx_println(cctx);
                flint_printf("nvars = %wd\n", mctx->nvars);
                flint_printf("f1 = "); gr_mpoly_print_pretty(f1, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f2 = "); gr_mpoly_print_pretty(f2, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("k1 = %wd\n", k1);
                flint_printf("k2 = %wd\n", k2);
                flint_printf("eq = %d %d %d %d\n", eq1, eq2, eq3, eq4);
                fflush(stdout);
                flint_abort();
            }

            if (n_randint(state, 2))
                status |= gr_mpoly_add(f1, f1, f2, mctx, cctx);
            else
                status |= gr_mpoly_mul(f1, f1, f2, mctx, cctx);

            eq1 = gr_mpoly_is_gen(f1, k1, mctx, cctx);
            eq2 = gr_mpoly_is_gen(f1, k2, mctx, cctx);
            eq3 = gr_mpoly_is_gen(f1, -1, mctx, cctx);

            if (status == GR_SUCCESS && (eq1 == T_TRUE || eq2 == T_TRUE || eq3 == T_TRUE))
            {
                flint_printf("FAIL\n");
                gr_ctx_println(cctx);
                flint_printf("nvars = %wd\n", mctx->nvars);
                flint_printf("f1 = "); gr_mpoly_print_pretty(f1, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f2 = "); gr_mpoly_print_pretty(f2, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("k1 = %wd\n", k1);
                flint_printf("k2 = %wd\n", k2);
                flint_printf("eq = %d %d %d\n", eq1, eq2, eq3);
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f1, mctx, cctx);
        gr_mpoly_clear(f2, mctx, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
