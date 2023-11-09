/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_div_monagan_pearce, state)
{
    int i, j, result;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k, l;
        mp_limb_t modulus;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);
        nmod_mpoly_init(l, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bits2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(l, state, len, exp_bits, ctx);

            nmod_mpoly_mul(h, f, g, ctx);

            nmod_mpoly_div_monagan_pearce(k, h, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            result = nmod_mpoly_equal(k, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_set(l, h, ctx);
            nmod_mpoly_div_monagan_pearce(l, l, g, ctx);
            nmod_mpoly_assert_canonical(l, ctx);
            result = nmod_mpoly_equal(l, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing dividend\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_div_monagan_pearce(g, h, g, ctx);
            nmod_mpoly_assert_canonical(g, ctx);
            result = nmod_mpoly_equal(g, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing divisor\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_clear(l, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check div matches divrem for random polys */
    for (i = 0; i < 40 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, q, r, k;
        mp_limb_t modulus;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong n;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(q, ctx);
        nmod_mpoly_init(r, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 50/n) + 1;
        exp_bound1 = n_randint(state, 50/n) + 1;
        exp_bound2 = n_randint(state, 50/n) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(q, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);

            nmod_mpoly_divrem_monagan_pearce(q, r, f, g, ctx);
            nmod_mpoly_assert_canonical(q, ctx);
            nmod_mpoly_assert_canonical(r, ctx);
            nmod_mpoly_remainder_strongtest(r, g, ctx);

            nmod_mpoly_mul_johnson(k, q, g, ctx);
            nmod_mpoly_add(k, k, r, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            result = nmod_mpoly_equal(f, k, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_div_monagan_pearce(k, f, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            result = nmod_mpoly_equal(k, q, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_set(k, f, ctx);
            nmod_mpoly_div_monagan_pearce(k, k, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            result = nmod_mpoly_equal(k, q, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem aliasing dividend\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            nmod_mpoly_div_monagan_pearce(g, f, g, ctx);
            nmod_mpoly_assert_canonical(g, ctx);
            result = nmod_mpoly_equal(g, q, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check div matches divrem aliasing divisor\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(q, ctx);
        nmod_mpoly_clear(r, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
