/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_mul_johnson, state)
{
    int i, j, result;
    slong tmul = 5;

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < 2 * tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k1, ctx);
        fq_nmod_mpoly_init(k2, ctx);
        fq_nmod_mpoly_init(t1, ctx);
        fq_nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        fq_nmod_mpoly_randtest_bits(k1, state, len, exp_bits, ctx);
        fq_nmod_mpoly_randtest_bits(k2, state, len, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            fq_nmod_mpoly_add(t1, g, h, ctx);
            fq_nmod_mpoly_assert_canonical(t1, ctx);
            fq_nmod_mpoly_mul_johnson(k1, f, t1, ctx);
            fq_nmod_mpoly_assert_canonical(k1, ctx);
            fq_nmod_mpoly_mul_johnson(t1, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(t1, ctx);
            fq_nmod_mpoly_mul_johnson(t2, f, h, ctx);
            fq_nmod_mpoly_assert_canonical(t2, ctx);
            fq_nmod_mpoly_add(k2, t1, t2, ctx);
            fq_nmod_mpoly_assert_canonical(k2, ctx);
            result = fq_nmod_mpoly_equal(k1, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*(g + h) = f*g + f*h\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(k1, ctx);
        fq_nmod_mpoly_clear(k2, ctx);
        fq_nmod_mpoly_clear(t1, ctx);
        fq_nmod_mpoly_clear(t2, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            fq_nmod_mpoly_mul_johnson(h, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_mul_johnson(f, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);
            result = fq_nmod_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            fq_nmod_mpoly_mul_johnson(h, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_mul_johnson(g, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);
            result = fq_nmod_mpoly_equal(h, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
