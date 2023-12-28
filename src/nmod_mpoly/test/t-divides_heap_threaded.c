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

#if defined(nmod_mpoly_divides_heap_threaded)
TEST_FUNCTION_START(nmod_mpoly_divides_heap_threaded, state)
{
    int i, j, result, result2, max_threads = 5, tmul = 30;
#ifdef _WIN32
    tmul = 1;
#endif

    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t p, f, g, h, h2;
        const char * vars[] = {"x","y","z","t","u"};
        int aff[] = {0, 1};

        nmod_mpoly_ctx_init(ctx, 5, ORD_DEGLEX, 179424691);
        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(h2, ctx);
        nmod_mpoly_set_str_pretty(g, "(1+x+y+2*z^2+3*t^3+5*u^5)^6", vars, ctx);
        nmod_mpoly_set_str_pretty(f, "(1+u+t+2*z^2+3*y^3+5*x^5)^6", vars, ctx);

        nmod_mpoly_mul(p, f, g, ctx);
        result = nmod_mpoly_divides_monagan_pearce(h, p, f, ctx);
        nmod_mpoly_assert_canonical(h, ctx);
        flint_set_num_threads(2);
        flint_set_thread_affinity(aff, 2);
        result2 = nmod_mpoly_divides_heap_threaded(h2, p, f, ctx);
        nmod_mpoly_assert_canonical(h2, ctx);
        flint_restore_thread_affinity();

        if (!result || !result2
                    || !nmod_mpoly_equal(h, g, ctx)
                    || !nmod_mpoly_equal(h2, g, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check simple example\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(h2, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            result = nmod_mpoly_divides_heap_threaded(k, h, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            result = result && nmod_mpoly_equal(f, k, ctx);

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
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f with divisor aliasing */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            result = nmod_mpoly_divides_heap_threaded(g, h, g, ctx);
            nmod_mpoly_assert_canonical(g, ctx);
            result = result && nmod_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f with divisor aliasing\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f with dividend aliasing */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            result = nmod_mpoly_divides_heap_threaded(h, h, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            result = result && nmod_mpoly_equal(f, h, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f with dividend aliasing\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check random polys don't divide */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, p, h1, h2;
        slong len1, len2, len3;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bound3;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(h1, ctx);
        nmod_mpoly_init(h2, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;
        len3 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;
        exp_bound3 = n_randint(state, 20) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(p, state, len3, exp_bound3, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul(f, f, g, ctx);
            nmod_mpoly_add(f, f, p, ctx);
            result = nmod_mpoly_divides_monagan_pearce(h1, f, g, ctx);
            nmod_mpoly_assert_canonical(h1, ctx);
            result2 = nmod_mpoly_divides_heap_threaded(h2, f, g, ctx);
            nmod_mpoly_assert_canonical(h2, ctx);

            if (result != result2 || !nmod_mpoly_equal(h1, h2, ctx))
            {
                flint_printf("Check random polys don't divide\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_clear(h1, ctx);
        nmod_mpoly_clear(h2, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check random polys don't divide alias dividend */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, p, h1;
        slong len1, len2, len3;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bound3;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(h1, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;
        len3 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;
        exp_bound3 = n_randint(state, 20) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(p, state, len3, exp_bound3, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul(f, f, g, ctx);
            nmod_mpoly_add(f, f, p, ctx);
            result = nmod_mpoly_divides_monagan_pearce(h1, f, g, ctx);
            nmod_mpoly_assert_canonical(h1, ctx);
            result2 = nmod_mpoly_divides_heap_threaded(f, f, g, ctx);
            nmod_mpoly_assert_canonical(f, ctx);

            if (result != result2 || !nmod_mpoly_equal(h1, f, ctx))
            {
                flint_printf("Check random polys don't divide alias dividend\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_clear(h1, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check random polys don't divide alias divisor */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, p, h1;
        slong len1, len2, len3;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bound3;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(p, ctx);
        nmod_mpoly_init(h1, ctx);

        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;
        len3 = n_randint(state, 10);

        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;
        exp_bound3 = n_randint(state, 20) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bound(p, state, len3, exp_bound3, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            nmod_mpoly_mul(f, f, g, ctx);
            nmod_mpoly_add(f, f, p, ctx);
            result = nmod_mpoly_divides_monagan_pearce(h1, f, g, ctx);
            nmod_mpoly_assert_canonical(h1, ctx);
            result2 = nmod_mpoly_divides_heap_threaded(g, f, g, ctx);
            nmod_mpoly_assert_canonical(g, ctx);

            if (result != result2 || !nmod_mpoly_equal(h1, g, ctx))
            {
                flint_printf("Check random polys don't divide alias divisor\n"
                                                   "i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(p, ctx);
        nmod_mpoly_clear(h1, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
#else
TEST_FUNCTION_START(nmod_mpoly_divides_heap_threaded, state)
{
    TEST_FUNCTION_END_SKIPPED(state);
}
#endif
