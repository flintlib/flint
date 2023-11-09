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

TEST_FUNCTION_START(nmod_mpoly_divides_dense, state)
{
    int i, j, result, ret;

    /* Check f*g/g = f */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        mp_limb_t modulus;
        slong n;

        modulus = n_randint(state, (i % 4 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 6, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 100/n/n;
        exp_bound = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
        }

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bounds(g, state, len2 + 1, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bounds(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            nmod_mpoly_mul(h, f, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            ret = nmod_mpoly_divides_dense(k, h, g, ctx);
            if (ret == -1)
                continue;
            nmod_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && nmod_mpoly_equal(k, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check divisibility of random polys */
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        mp_limb_t modulus;
        slong n;

        modulus = n_randint(state, (i % 4 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 6, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 20/n;
        exp_bound = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
        }

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bounds(g, state, len2 + 1, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bounds(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            ret = nmod_mpoly_divides_dense(h, f, g, ctx);
            if (ret == -1 || ret == 0)
                continue;
            nmod_mpoly_assert_canonical(h, ctx);
            nmod_mpoly_mul(k, h, g, ctx);
            nmod_mpoly_assert_canonical(k, ctx);
            result = nmod_mpoly_equal(k, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility of random polys\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f aliasing first argument */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        mp_limb_t modulus;
        slong n;

        modulus = n_randint(state, (i % 4 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 6, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 100/n/n;
        exp_bound = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
        }

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bounds(g, state, len2 + 1, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bounds(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            nmod_mpoly_set(h, f, ctx);
            nmod_mpoly_mul(h, h, g, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            ret = nmod_mpoly_divides_dense(k, h, g, ctx);
            if (ret == -1)
                continue;
            nmod_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && nmod_mpoly_equal(k, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f aliasing second argument */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        mp_limb_t modulus;
        slong n;

        modulus = n_randint(state, (i % 4 == 0) ? 4: FLINT_BITS - 1) + 1;
        modulus = n_randbits(state, modulus);
        modulus = n_nextprime(modulus, 1);

        nmod_mpoly_ctx_init_rand(ctx, state, 6, modulus);

        nmod_mpoly_init(f, ctx);
        nmod_mpoly_init(g, ctx);
        nmod_mpoly_init(h, ctx);
        nmod_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 100/n/n;
        exp_bound = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
        }

        for (j = 0; j < 4; j++)
        {
            nmod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            nmod_mpoly_randtest_bounds(g, state, len2 + 1, exp_bound2, ctx);
            if (nmod_mpoly_is_zero(g, ctx))
                nmod_mpoly_one(g, ctx);
            nmod_mpoly_randtest_bounds(h, state, len, exp_bound, ctx);
            nmod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            nmod_mpoly_set(h, g, ctx);
            nmod_mpoly_mul(h, f, h, ctx);
            nmod_mpoly_assert_canonical(h, ctx);
            ret = nmod_mpoly_divides_dense(k, h, g, ctx);
            if (ret == -1)
                continue;
            nmod_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && nmod_mpoly_equal(k, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing second argument\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_clear(g, ctx);
        nmod_mpoly_clear(h, ctx);
        nmod_mpoly_clear(k, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
