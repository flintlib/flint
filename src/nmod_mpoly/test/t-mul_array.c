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

TEST_FUNCTION_START(nmod_mpoly_mul_array, state)
{
    int i, j, result;

    /* Check mul_array matches mul_johnson */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong n, max_bound;
        mp_limb_t modulus;

        modulus = n_randint(state, -UWORD(2)) + UWORD(2);
        nmod_mpoly_ctx_init_rand(ctx, state, 5, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        max_bound = ctx->minfo->ord == ORD_LEX ? 200 : 100;
        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = max_bound/n/n;
        exp_bound =  n_randint(state, max_bound) + 1;
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            result = nmod_mpoly_mul_array(k, f, g, ctx);
            if (!result)
            {
                continue;
            }
            nmod_mpoly_assert_canonical(k, ctx);
            result = nmod_mpoly_equal(h, k, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check mul_array matches mul_johnson\ni = %wd, j = %wd\n", i, j);
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

    /* Check aliasing first argument */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong n, max_bound;
        mp_limb_t modulus;

        modulus = n_randint(state, -UWORD(2)) + UWORD(2);
        nmod_mpoly_ctx_init_rand(ctx, state, 5, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 200/n/n;
        exp_bound =  n_randint(state, max_bound) + 1;
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            result = nmod_mpoly_mul_array(f, f, g, ctx);
            if (!result)
                continue;

            nmod_mpoly_assert_canonical(f, ctx);
            result = nmod_mpoly_equal(h, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing first argument\ni = %wd, j = %wd\n", i, j);
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
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong n, max_bound;
        mp_limb_t modulus;

        modulus = n_randint(state, -UWORD(2)) + UWORD(2);
        nmod_mpoly_ctx_init_rand(ctx, state, 5, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 200/n/n;
        exp_bound =  n_randint(state, max_bound) + 1;
        exp_bound1 = n_randint(state, max_bound) + 1;
        exp_bound2 = n_randint(state, max_bound) + 1;

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);

            nmod_mpoly_mul_johnson(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            result = nmod_mpoly_mul_array(f, g, f, ctx);
            if (!result)
                continue;

            nmod_mpoly_assert_canonical(f, ctx);
            result = nmod_mpoly_equal(h, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second argument\ni = %wd, j = %wd\n", i, j);
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
