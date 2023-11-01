/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_mul, state)
{
    int i, j, result, max_threads = 5;
    slong tmul = 10;
#ifdef _WIN32
    tmul = 2;
#endif

    /* check fixed cases */
    for (i = 0; i < 1 + tmul; i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        const char * vars[] = {"x", "y", "z", "t"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_DEGLEX);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        /* should trigger dense */
        fmpz_mpoly_set_str_pretty(f, "((1-x)*(1+y)*(1+z))^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "((1+x)*(1-y)*(1-z))^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(k, "((1-x^2)*(1-y^2)*(1-z^2))^20", vars, ctx);
        fmpz_mpoly_mul(h, f, g, ctx);
        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check fixed case 1\n");
            fflush(stdout);
            flint_abort();
        }

        /* should trigger array */
        fmpz_mpoly_set_str_pretty(f, "(1+x+y+z+t)^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "(1-x-y-z-t)^20", vars, ctx);
        fmpz_mpoly_set_str_pretty(k, "((1+x+y+z+t)*(1-x-y-z-t))^20", vars, ctx);
        fmpz_mpoly_mul(h, f, g, ctx);
        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check fixed case 2\n");
            fflush(stdout);
            flint_abort();
        }

        /* should trigger heap */
        fmpz_mpoly_set_str_pretty(f, "(1+x^10)^50", vars, ctx);
        fmpz_mpoly_set_str_pretty(g, "(1+y^10)^50", vars, ctx);
        fmpz_mpoly_set_str_pretty(k, "((1+x^10)*(1+y^10))^50", vars, ctx);
        fmpz_mpoly_mul(h, f, g, ctx);
        if (!fmpz_mpoly_equal(h, k, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check fixed case 3\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);

        flint_set_num_threads(n_randint(state, max_threads) + 1);
    }

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 20);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k1, ctx);
        fmpz_mpoly_init(k2, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            exp_bits = n_randint(state, 100) + 2;
            exp_bits1 = n_randint(state, 100) + 2;
            exp_bits2 = n_randint(state, 100) + 2;

            fmpz_mpoly_randtest_bits(k1, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k2, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpz_mpoly_randtest_bits(h, state, len2, coeff_bits, exp_bits2, ctx);

            fmpz_mpoly_add(t1, g, h, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            fmpz_mpoly_mul(k1, f, t1, ctx);
            fmpz_mpoly_assert_canonical(k1, ctx);
            fmpz_mpoly_mul(t1, f, g, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            flint_set_num_threads(n_randint(state, max_threads) + 1);
            fmpz_mpoly_mul(t2, f, h, ctx);
            fmpz_mpoly_assert_canonical(t2, ctx);
            fmpz_mpoly_add(k2, t1, t2, ctx);
            fmpz_mpoly_assert_canonical(k2, ctx);
            result = fmpz_mpoly_equal(k1, k2, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*(g + h) = f*g + f*h\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k1, ctx);
        fmpz_mpoly_clear(k2, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1, len2, n;
        flint_bitcnt_t coeff_bits, exp_bound, exp_bound1, exp_bound2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        exp_bound = 3 + n_randint(state, 1 + 100/n/n);
        exp_bound1 = 3 + n_randint(state, 1 + 100/n/n);
        exp_bound2 = 3 + n_randint(state, 1 + 100/n/n);

        len = exp_bound + 1;
        len1 = exp_bound1 + 1;
        len2 = exp_bound2 + 1;
        for (j = n_randint(state, n) + 2; j >= 0; j--)
        {
            len *= exp_bound + 1;
            len1 *= exp_bound1 + 1;
            len2 *= exp_bound2 + 1;
            len = FLINT_MIN(len, WORD(1000));
            len1 = FLINT_MIN(len, WORD(1000));
            len2 = FLINT_MIN(len, WORD(1000));
        }

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            flint_set_num_threads(n_randint(state, max_threads) + 1);
            fmpz_mpoly_mul(f, f, g, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            result = fmpz_mpoly_equal(h, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h;
        slong len, len1, len2, n;
        flint_bitcnt_t coeff_bits, exp_bound, exp_bound1, exp_bound2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 4);
        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);

        exp_bound = 3 + n_randint(state, 1 + 100/n/n);
        exp_bound1 = 3 + n_randint(state, 1 + 100/n/n);
        exp_bound2 = 3 + n_randint(state, 1 + 100/n/n);

        len = exp_bound + 1;
        len1 = exp_bound1 + 1;
        len2 = exp_bound2 + 1;
        for (j = n_randint(state, n) + 2; j >= 0; j--)
        {
            len *= exp_bound + 1;
            len1 *= exp_bound1 + 1;
            len2 *= exp_bound2 + 1;
            len = FLINT_MIN(len, WORD(1000));
            len1 = FLINT_MIN(len, WORD(1000));
            len2 = FLINT_MIN(len, WORD(1000));
        }

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 2; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits, exp_bound2, ctx);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            flint_set_num_threads(n_randint(state, max_threads) + 1);
            fmpz_mpoly_mul(g, f, g, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            result = fmpz_mpoly_equal(h, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
