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

TEST_FUNCTION_START(nmod_mpoly_divrem_monagan_pearce, state)
{
    int i, j, result;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k, r;
        mp_limb_t modulus;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);
        nmod_mpoly_init(r, ctx);

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
            nmod_mpoly_randtest_bits(r, state, len, exp_bits, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);

            nmod_mpoly_divrem_monagan_pearce(k, r, h, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            nmod_mpoly_assert_canonical(r, ctx);

            result = nmod_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_clear(r, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f = g*q + r for random polys */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k, r;
        mp_limb_t modulus;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);
        nmod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound2 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);

            nmod_mpoly_divrem_monagan_pearce(h, r, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_assert_canonical(r, ctx);
            nmod_mpoly_remainder_strongtest(r, g, ctx);

            nmod_mpoly_mul_johnson(k, h, g, ctx);
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
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_clear(r, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of quotient with first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, r1, r2;
        mp_limb_t modulus;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(r1, ctx);
        nmod_mpoly_init(r2, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound2 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(r1, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(r2, state, len, exp_bound, ctx);

            nmod_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_assert_canonical(r1, ctx);
            nmod_mpoly_remainder_strongtest(r1, g, ctx);
            nmod_mpoly_divrem_monagan_pearce(f, r2, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_assert_canonical(r2, ctx);
            nmod_mpoly_remainder_strongtest(r2, g, ctx);

            result = nmod_mpoly_equal(h, f, ctx) && nmod_mpoly_equal(r1, r2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of quotient with first argument\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(r1, ctx);
        nmod_mpoly_clear(r2, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of quotient with second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, r1, r2;
        mp_limb_t modulus;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(r1, ctx);
        nmod_mpoly_init(r2, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound2 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(r1, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(r2, state, len, exp_bound, ctx);

            nmod_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_assert_canonical(r1, ctx);
            nmod_mpoly_remainder_strongtest(r1, g, ctx);
            nmod_mpoly_divrem_monagan_pearce(g, r2, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_assert_canonical(r2, ctx);

            result = nmod_mpoly_equal(h, g, ctx) && nmod_mpoly_equal(r1, r2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of quotient with second argument\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(r1, ctx);
        nmod_mpoly_clear(r2, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of remainder with first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k, r1;
        mp_limb_t modulus;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);
        nmod_mpoly_init(r1, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound2 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(r1, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);

            nmod_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_assert_canonical(r1, ctx);
            nmod_mpoly_remainder_strongtest(r1, g, ctx);

            nmod_mpoly_divrem_monagan_pearce(k, f, f, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            nmod_mpoly_assert_canonical(f, ctx);
            nmod_mpoly_remainder_strongtest(f, g, ctx);

            result = nmod_mpoly_equal(h, k, ctx) && nmod_mpoly_equal(r1, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of remainder with first argument\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_clear(r1, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of remainder with second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k, r1;
        mp_limb_t modulus;
        slong nvars, len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);
        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);
        nvars = ctx->minfo->nvars;

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);
        nmod_mpoly_init(r1, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound1 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;
        exp_bound2 = n_randint(state, 50/FLINT_MAX(WORD(1), nvars)) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(r1, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);

            nmod_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_assert_canonical(r1, ctx);
            nmod_mpoly_remainder_strongtest(r1, g, ctx);

            nmod_mpoly_divrem_monagan_pearce(k, g, f, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            nmod_mpoly_assert_canonical(g, ctx);

            result = nmod_mpoly_equal(h, k, ctx) && nmod_mpoly_equal(r1, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of remainder with second argument\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_clear(r1, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
