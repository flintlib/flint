/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_mul, state)
{
    int i, j, result, max_threads = 5;
    int tmul = 10;
#ifdef _WIN32
    tmul = 1;
#endif

    /* Check f*(g + h) = f*g + f*h with bit bound */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, SMALL_FMPZ_BITCOUNT_MAX) + 2;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx, state, 5, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k1, ctx);
        nmod_mpoly_init(k2, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 100) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        nmod_mpoly_randtest_bits(k1, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(k2, state, len, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len2, exp_bits, ctx);

            nmod_mpoly_add(t1, g, h, ctx);
            nmod_mpoly_assert_canonical(t1, ctx);
            nmod_mpoly_mul(k1, f, t1, ctx);
            nmod_mpoly_assert_canonical(k1, ctx);
            nmod_mpoly_mul(t1, f, g, ctx);
            nmod_mpoly_assert_canonical(t1, ctx);
            nmod_mpoly_mul(t2, f, h, ctx);
            nmod_mpoly_assert_canonical(t2, ctx);
            nmod_mpoly_add(k2, t1, t2, ctx);
            nmod_mpoly_assert_canonical(k2, ctx);
            result = nmod_mpoly_equal(k1, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*(g + h) = f*g + f*h with bit bound\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k1, ctx);
        nmod_mpoly_clear(k2, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t exp_bound, exp_bound1, exp_bound2;
        mp_limb_t modulus;
        slong n;

        modulus = n_randint(state, SMALL_FMPZ_BITCOUNT_MAX) + 2;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx, state, 4, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
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

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_mul(f, f, g, ctx);
            nmod_mpoly_assert_canonical(f, ctx);
            result = nmod_mpoly_equal(h, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t exp_bound, exp_bound1, exp_bound2;
        mp_limb_t modulus;
        slong n;

        modulus = n_randint(state, SMALL_FMPZ_BITCOUNT_MAX) + 2;
        modulus = n_randbits(state, modulus);
        nmod_mpoly_ctx_init_rand(ctx, state, 4, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
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

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_mul(g, f, g, ctx);
            nmod_mpoly_assert_canonical(g, ctx);
            result = nmod_mpoly_equal(h, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
