/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpoly.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_add_sub, state)
{
    slong i, j;

    for (i = 0; i < 100; i++)
    {
        gr_ctx_t cctx;
        gr_mpoly_ctx_t ctx;
        gr_mpoly_t f, g, h, k, k1, k2;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        int status;

        gr_ctx_init_random(cctx, state);
        gr_mpoly_ctx_init_rand(ctx, state, cctx, 10);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(g, ctx);
        gr_mpoly_init(h, ctx);
        gr_mpoly_init(k, ctx);
        gr_mpoly_init(k1, ctx);
        gr_mpoly_init(k2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        /* Check (f + g) - g = f */
        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            status |= gr_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

            status |= gr_mpoly_add(h, g, f, ctx);
            status |= gr_mpoly_sub(k, h, g, ctx);

            if (status == GR_SUCCESS)
            {
                gr_mpoly_assert_canonical(h, ctx);
                gr_mpoly_assert_canonical(k, ctx);
            }

            if (status == GR_SUCCESS && gr_mpoly_equal(f, k, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Check (f + g) - g = f\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                flint_printf("f + g = "); gr_mpoly_print_pretty(h, ctx); flint_printf("\n");
                flint_printf("(f + g) - g = "); gr_mpoly_print_pretty(k, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check f + g = g + f */
        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            status |= gr_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

            status |= gr_mpoly_add(h, f, g, ctx);
            status |= gr_mpoly_add(k, g, f, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(h, k, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Check (f + g) = (g + f)\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                flint_printf("f + g = "); gr_mpoly_print_pretty(h, ctx); flint_printf("\n");
                flint_printf("g + f = "); gr_mpoly_print_pretty(k, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check f - g = -g + f */
        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            status |= gr_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

            status |= gr_mpoly_sub(h, f, g, ctx);
            status |= gr_mpoly_neg(k, g, ctx);
            status |= gr_mpoly_add(k, k, f, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(h, k, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Check f - g = -g + f\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                flint_printf("f - g = "); gr_mpoly_print_pretty(h, ctx); flint_printf("\n");
                flint_printf("-g + f = "); gr_mpoly_print_pretty(k, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check f + (g + h) = (f + g) + h */
        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            status |= gr_mpoly_randtest_bits(k1, state, len, exp_bits, ctx);
            status |= gr_mpoly_randtest_bits(k2, state, len, exp_bits, ctx);

            status |= gr_mpoly_add(k1, g, h, ctx);
            status |= gr_mpoly_add(k1, f, k1, ctx);
            status |= gr_mpoly_add(k2, f, g, ctx);
            status |= gr_mpoly_add(k2, k2, h, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(k1, k2, ctx) == T_FALSE)
            {
                flint_printf("FAIL Check f + (g + h) = (f + g) + h\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                flint_printf("h = "); gr_mpoly_print_pretty(h, ctx); flint_printf("\n");
                flint_printf("f + (g + h) = "); gr_mpoly_print_pretty(k1, ctx); flint_printf("\n");
                flint_printf("(f + g) + h = "); gr_mpoly_print_pretty(k2, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check f - (g + h) = (f - g) - h */
        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            status |= gr_mpoly_randtest_bits(k1, state, len, exp_bits, ctx);
            status |= gr_mpoly_randtest_bits(k2, state, len, exp_bits, ctx);

            status |= gr_mpoly_add(k1, g, h, ctx);
            status |= gr_mpoly_sub(k1, f, k1, ctx);
            status |= gr_mpoly_sub(k2, f, g, ctx);
            status |= gr_mpoly_sub(k2, k2, h, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(k1, k2, ctx) == T_FALSE)
            {
                flint_printf("FAIL Check f + (g + h) = (f + g) + h\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(ctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, ctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, ctx); flint_printf("\n");
                flint_printf("h = "); gr_mpoly_print_pretty(h, ctx); flint_printf("\n");
                flint_printf("f + (g + h) = "); gr_mpoly_print_pretty(k1, ctx); flint_printf("\n");
                flint_printf("(f + g) + h = "); gr_mpoly_print_pretty(k2, ctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        /* Check aliasing */
        for (j = 0; j < 1; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            status |= gr_mpoly_set(h, f, ctx);

            status |= gr_mpoly_add(f, f, g, ctx);
            gr_mpoly_assert_canonical(f, ctx);
            status |= gr_mpoly_sub(f, f, g, ctx);
            gr_mpoly_assert_canonical(f, ctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(f, h, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Check aliasing first arg\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            if ((j % 2) == 0)
            {
                status |= gr_mpoly_add(h, g, f, ctx);
                gr_mpoly_assert_canonical(h, ctx);
                status |= gr_mpoly_add(f, g, f, ctx);
                gr_mpoly_assert_canonical(f, ctx);
            }
            else
            {
                status |= gr_mpoly_sub(h, g, f, ctx);
                gr_mpoly_assert_canonical(h, ctx);
                status |= gr_mpoly_sub(f, g, f, ctx);
                gr_mpoly_assert_canonical(f, ctx);
            }

            if (status == GR_SUCCESS && gr_mpoly_equal(f, h, ctx) == T_FALSE)
            {
                flint_printf("FAIL: Check aliasing second arg\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, ctx);
        gr_mpoly_clear(g, ctx);
        gr_mpoly_clear(h, ctx);
        gr_mpoly_clear(k, ctx);
        gr_mpoly_clear(k1, ctx);
        gr_mpoly_clear(k2, ctx);

        gr_mpoly_ctx_clear(ctx);
        gr_ctx_clear(cctx);
    }

    TEST_FUNCTION_END(state);
}
