/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "nmod_mpoly.h"

int
main(void)
{
    int i, j, result;
    FLINT_TEST_INIT(state);

    flint_printf("mul....");
    fflush(stdout);

    /* Check f*(g + h) = f*g + f*h with bit bound */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 5, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k1, ctx);
        nmod_mpoly_init(k2, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 200);
        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        exp_bits = n_randint(state, 20) + 1;
        exp_bits1 = n_randint(state, 20) + 1;
        exp_bits2 = n_randint(state, 20) + 1;

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

    /* Check f*(g + h) = f*g + f*h with size bound */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        slong max_bound, exp_bound, exp_bound1, exp_bound2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 5, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k1, ctx);
        nmod_mpoly_init(k2, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 200);
        len1 = n_randint(state, 200);
        len2 = n_randint(state, 200);

        max_bound = 1 + 1000/ctx->minfo->nvars;
        exp_bound = n_randint(state, max_bound) + 1;
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        nmod_mpoly_randtest_bound(k1, state, len, exp_bound, ctx);
        nmod_mpoly_randtest_bound(k2, state, len, exp_bound, ctx);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest_bound(h, state, len2, exp_bound, ctx);

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
                flint_printf("Check f*(g + h) = f*g + f*h with size bound\ni = %wd, j = %wd\n", i ,j);
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
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        slong len, len1, len2;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 100) + 1;
        exp_bits1 = n_randint(state, 100) + 1;
        exp_bits2 = n_randint(state, 100) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            nmod_mpoly_mul(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_mul(f, f, g, ctx);
            nmod_mpoly_assert_canonical(f, ctx);
            result = nmod_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        slong len, len1, len2;
        mp_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        mp_limb_t modulus;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 10, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 100) + 2;
        exp_bits1 = n_randint(state, 100) + 2;
        exp_bits2 = n_randint(state, 100) + 2;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            nmod_mpoly_mul(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_mul(g, f, g, ctx);
            nmod_mpoly_assert_canonical(g, ctx);
            result = nmod_mpoly_equal(h, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

