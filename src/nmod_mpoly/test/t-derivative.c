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

TEST_FUNCTION_START(nmod_mpoly_derivative, state)
{
    int i, j, result;
    slong tmul = 5;

    /* Check d(f*g) = df*g + f*dg */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, fp, gp, hp, t1, t2;
        mp_limb_t modulus;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        slong idx;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        if (ctx->minfo->nvars < 1)
        {
            nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(fp, ctx);
        nmod_mpoly_init(gp, ctx);
        nmod_mpoly_init(hp, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        nmod_mpoly_randtest_bits(hp, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(fp, state, len, exp_bits1, ctx);
        nmod_mpoly_randtest_bits(gp, state, len, exp_bits2, ctx);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bits2, ctx);

            idx = n_randint(state, ctx->minfo->nvars);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);

            nmod_mpoly_derivative(hp, h, idx, ctx);
            nmod_mpoly_assert_canonical(hp, ctx);

            nmod_mpoly_derivative(fp, f, idx, ctx);
            nmod_mpoly_assert_canonical(fp, ctx);
            nmod_mpoly_derivative(gp, g, idx, ctx);
            nmod_mpoly_assert_canonical(gp, ctx);

            nmod_mpoly_mul_johnson(t1, f, gp, ctx);
            nmod_mpoly_assert_canonical(t1, ctx);
            nmod_mpoly_mul_johnson(t2, g, fp, ctx);
            nmod_mpoly_assert_canonical(t2, ctx);
            nmod_mpoly_add(t1, t1, t2, ctx);
            nmod_mpoly_assert_canonical(t1, ctx);

            result = nmod_mpoly_equal(hp, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(fp, ctx);
        nmod_mpoly_clear(gp, ctx);
        nmod_mpoly_clear(hp, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check d(f*g) = df*g + f*dg with aliasing */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, fp, gp, t1, t2;
        mp_limb_t modulus;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        slong idx;

        modulus = n_randint(state, FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);

        nmod_mpoly_ctx_init_rand(ctx, state, 20, modulus);
        if (ctx->minfo->nvars < 1)
        {
            nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(fp, ctx);
        nmod_mpoly_init(gp, ctx);
        nmod_mpoly_init(t1, ctx);
        nmod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(fp, state, len, exp_bits, ctx);
        nmod_mpoly_randtest_bits(gp, state, len, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);

            idx = n_randint(state, ctx->minfo->nvars);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);

            nmod_mpoly_derivative(h, h, idx, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_set(fp, f, ctx);
            nmod_mpoly_derivative(fp, fp, idx, ctx);
            nmod_mpoly_assert_canonical(fp, ctx);
            nmod_mpoly_set(gp, g, ctx);
            nmod_mpoly_derivative(gp, gp, idx, ctx);
            nmod_mpoly_assert_canonical(gp, ctx);

            nmod_mpoly_mul_johnson(t1, f, gp, ctx);
            nmod_mpoly_assert_canonical(t1, ctx);
            nmod_mpoly_mul_johnson(t2, g, fp, ctx);
            nmod_mpoly_assert_canonical(t2, ctx);
            nmod_mpoly_add(t1, t1, t2, ctx);
            nmod_mpoly_assert_canonical(t1, ctx);

            result = nmod_mpoly_equal(h, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg with aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(fp, ctx);
        nmod_mpoly_clear(gp, ctx);
        nmod_mpoly_clear(t1, ctx);
        nmod_mpoly_clear(t2, ctx);

        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
