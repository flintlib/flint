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

TEST_FUNCTION_START(fq_nmod_mpoly_div_monagan_pearce, state)
{
    int i, j, result;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k, l;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);
        fq_nmod_mpoly_init(l, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bits2, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(l, state, len, exp_bits, ctx);

            fq_nmod_mpoly_mul(h, f, g, ctx);

            fq_nmod_mpoly_div_monagan_pearce(k, h, g, ctx);
            fq_nmod_mpoly_assert_canonical(k, ctx);
            result = fq_nmod_mpoly_equal(k, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_set(l, h, ctx);
            fq_nmod_mpoly_div_monagan_pearce(l, l, g, ctx);
            fq_nmod_mpoly_assert_canonical(l, ctx);
            result = fq_nmod_mpoly_equal(l, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing dividend\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_div_monagan_pearce(g, h, g, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);
            result = fq_nmod_mpoly_equal(g, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing divisor\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(k, ctx);
        fq_nmod_mpoly_clear(l, ctx);

        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check div matches divrem for random polys */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, q, r, k;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong n;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(q, ctx);
        fq_nmod_mpoly_init(r, ctx);
        fq_nmod_mpoly_init(k, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 15);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 100/n) + 1;
        exp_bound1 = n_randint(state, 100/n) + 1;
        exp_bound2 = n_randint(state, 100/n) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bound(q, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(q, r, f, g, ctx);

            fq_nmod_mpoly_assert_canonical(q, ctx);
            fq_nmod_mpoly_assert_canonical(r, ctx);
            fq_nmod_mpoly_remainder_strongtest(r, g, ctx);

            fq_nmod_mpoly_mul_johnson(k, q, g, ctx);
            fq_nmod_mpoly_add(k, k, r, ctx);
            fq_nmod_mpoly_assert_canonical(k, ctx);
            result = fq_nmod_mpoly_equal(f, k, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_div_monagan_pearce(k, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(k, ctx);
            result = fq_nmod_mpoly_equal(k, q, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_set(k, f, ctx);
            fq_nmod_mpoly_div_monagan_pearce(k, k, g, ctx);
            fq_nmod_mpoly_assert_canonical(k, ctx);
            result = fq_nmod_mpoly_equal(k, q, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem aliasing dividend\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_div_monagan_pearce(g, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);
            result = fq_nmod_mpoly_equal(g, q, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem aliasing divisor\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(q, ctx);
        fq_nmod_mpoly_clear(r, ctx);
        fq_nmod_mpoly_clear(k, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
