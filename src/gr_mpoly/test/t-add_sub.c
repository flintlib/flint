/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "gr_mpoly.h"

TEST_FUNCTION_START(gr_mpoly_add_sub, state)
{
    slong i, j;

    /* Check (f + g) - g = f */
    for (i = 0; i < 100; i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        int status;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 20);

        gr_mpoly_init(f, mctx, cctx);
        gr_mpoly_init(g, mctx, cctx);
        gr_mpoly_init(h, mctx, cctx);
        gr_mpoly_init(k, mctx, cctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, mctx, cctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, mctx, cctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k, state, len, exp_bits, mctx, cctx);

            status |= gr_mpoly_add(h, g, f, mctx, cctx);
            status |= gr_mpoly_sub(k, h, g, mctx, cctx);

            if (status == GR_SUCCESS)
            {
                gr_mpoly_assert_canonical(h, mctx, cctx);
                gr_mpoly_assert_canonical(k, mctx, cctx);
            }

            if (status == GR_SUCCESS && gr_mpoly_equal(f, k, mctx, cctx) == T_FALSE)
            {
                flint_printf("FAIL: Check (f + g) - g = f\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f + g = "); gr_mpoly_print_pretty(h, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("(f + g) - g = "); gr_mpoly_print_pretty(k, NULL, mctx, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, mctx, cctx);
        gr_mpoly_clear(g, mctx, cctx);
        gr_mpoly_clear(h, mctx, cctx);
        gr_mpoly_clear(k, mctx, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

    /* Check f + g = g + f */
    for (i = 0; i < 10; i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        int status;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 20);

        gr_mpoly_init(f, mctx, cctx);
        gr_mpoly_init(g, mctx, cctx);
        gr_mpoly_init(h, mctx, cctx);
        gr_mpoly_init(k, mctx, cctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, mctx, cctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, mctx, cctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k, state, len, exp_bits, mctx, cctx);

            status |= gr_mpoly_add(h, f, g, mctx, cctx);
            status |= gr_mpoly_add(k, g, f, mctx, cctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(h, k, mctx, cctx) == T_FALSE)
            {
                flint_printf("FAIL: Check (f + g) = (g + f)\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f + g = "); gr_mpoly_print_pretty(h, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("g + f = "); gr_mpoly_print_pretty(k, NULL, mctx, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, mctx, cctx);
        gr_mpoly_clear(g, mctx, cctx);
        gr_mpoly_clear(h, mctx, cctx);
        gr_mpoly_clear(k, mctx, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

    /* Check f - g = -g + f */
    for (i = 0; i < 10; i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        int status;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 20);

        gr_mpoly_init(f, mctx, cctx);
        gr_mpoly_init(g, mctx, cctx);
        gr_mpoly_init(h, mctx, cctx);
        gr_mpoly_init(k, mctx, cctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, mctx, cctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, mctx, cctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k, state, len, exp_bits, mctx, cctx);

            status |= gr_mpoly_sub(h, f, g, mctx, cctx);
            status |= gr_mpoly_neg(k, g, mctx, cctx);
            status |= gr_mpoly_add(k, k, f, mctx, cctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(h, k, mctx, cctx) == T_FALSE)
            {
                flint_printf("FAIL: Check f - g = -g + f\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f - g = "); gr_mpoly_print_pretty(h, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("-g + f = "); gr_mpoly_print_pretty(k, NULL, mctx, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, mctx, cctx);
        gr_mpoly_clear(g, mctx, cctx);
        gr_mpoly_clear(h, mctx, cctx);
        gr_mpoly_clear(k, mctx, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

    /* Check f + (g + h) = (f + g) + h */
    for (i = 0; i < 10; i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f, g, h, k1, k2;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        int status;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 20);

        gr_mpoly_init(f, mctx, cctx);
        gr_mpoly_init(g, mctx, cctx);
        gr_mpoly_init(h, mctx, cctx);
        gr_mpoly_init(k1, mctx, cctx);
        gr_mpoly_init(k2, mctx, cctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, mctx, cctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, mctx, cctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k1, state, len, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k2, state, len, exp_bits, mctx, cctx);

            status |= gr_mpoly_add(k1, g, h, mctx, cctx);
            status |= gr_mpoly_add(k1, f, k1, mctx, cctx);
            status |= gr_mpoly_add(k2, f, g, mctx, cctx);
            status |= gr_mpoly_add(k2, k2, h, mctx, cctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(k1, k2, mctx, cctx) == T_FALSE)
            {
                flint_printf("FAIL Check f + (g + h) = (f + g) + h\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("h = "); gr_mpoly_print_pretty(h, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f + (g + h) = "); gr_mpoly_print_pretty(k1, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("(f + g) + h = "); gr_mpoly_print_pretty(k2, NULL, mctx, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, mctx, cctx);
        gr_mpoly_clear(g, mctx, cctx);
        gr_mpoly_clear(h, mctx, cctx);
        gr_mpoly_clear(k1, mctx, cctx);
        gr_mpoly_clear(k2, mctx, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

    /* Check f - (g + h) = (f - g) - h */
    for (i = 0; i < 100; i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f, g, h, k1, k2;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        int status;

        gr_ctx_init_random(cctx, state);
        mpoly_ctx_init_rand(mctx, state, 20);

        gr_mpoly_init(f, mctx, cctx);
        gr_mpoly_init(g, mctx, cctx);
        gr_mpoly_init(h, mctx, cctx);
        gr_mpoly_init(k1, mctx, cctx);
        gr_mpoly_init(k2, mctx, cctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            status = GR_SUCCESS;

            status |= gr_mpoly_randtest_bits(f, state, len1, exp_bits1, mctx, cctx);
            status |= gr_mpoly_randtest_bits(g, state, len2, exp_bits2, mctx, cctx);
            status |= gr_mpoly_randtest_bits(h, state, len, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k1, state, len, exp_bits, mctx, cctx);
            status |= gr_mpoly_randtest_bits(k2, state, len, exp_bits, mctx, cctx);

            status |= gr_mpoly_add(k1, g, h, mctx, cctx);
            status |= gr_mpoly_sub(k1, f, k1, mctx, cctx);
            status |= gr_mpoly_sub(k2, f, g, mctx, cctx);
            status |= gr_mpoly_sub(k2, k2, h, mctx, cctx);

            if (status == GR_SUCCESS && gr_mpoly_equal(k1, k2, mctx, cctx) == T_FALSE)
            {
                flint_printf("FAIL Check f + (g + h) = (f + g) + h\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                gr_ctx_println(cctx);
                flint_printf("f = "); gr_mpoly_print_pretty(f, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("g = "); gr_mpoly_print_pretty(g, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("h = "); gr_mpoly_print_pretty(h, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("f + (g + h) = "); gr_mpoly_print_pretty(k1, NULL, mctx, cctx); flint_printf("\n");
                flint_printf("(f + g) + h = "); gr_mpoly_print_pretty(k2, NULL, mctx, cctx); flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        gr_mpoly_clear(f, mctx, cctx);
        gr_mpoly_clear(g, mctx, cctx);
        gr_mpoly_clear(h, mctx, cctx);
        gr_mpoly_clear(k1, mctx, cctx);
        gr_mpoly_clear(k2, mctx, cctx);

        mpoly_ctx_clear(mctx);
        gr_ctx_clear(cctx);
    }

#if 0
    /* Check aliasing */
    for (i = 0; i < 100; i++)
    {
        gr_ctx_t cctx;
        mpoly_ctx_t mctx;
        gr_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        gr_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);

        gr_mpoly_init(f, ctx);
        gr_mpoly_init(g, ctx);
        gr_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            gr_mpoly_set(h, f, ctx);

            gr_mpoly_add(f, f, g, ctx);
            gr_mpoly_assert_canonical(f, ctx);
            gr_mpoly_sub(f, f, g, ctx);
            gr_mpoly_assert_canonical(f, ctx);
            if (!gr_mpoly_equal(f, h, ctx))
            {
                flint_printf("FAIL: Check aliasing first arg\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            gr_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            gr_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            gr_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            if ((j % 2) == 0)
            {
                gr_mpoly_add(h, g, f, ctx);
                gr_mpoly_assert_canonical(h, ctx);
                gr_mpoly_add(f, g, f, ctx);
                gr_mpoly_assert_canonical(f, ctx);
            }
            else
            {
                gr_mpoly_sub(h, g, f, ctx);
                gr_mpoly_assert_canonical(h, ctx);
                gr_mpoly_sub(f, g, f, ctx);
                gr_mpoly_assert_canonical(f, ctx);
            }

            if (!gr_mpoly_equal(f, h, ctx))
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

        gr_mpoly_ctx_clear(ctx);
    }
#endif

    TEST_FUNCTION_END(state);
}
