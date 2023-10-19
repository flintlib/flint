/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_divides, state)
{
    slong i, j, result, ret, max_threads = 5, tmul = 25;

    /* Check f*g/g = f sparse */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k, hsave, gsave;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 6, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);
        fmpz_mod_mpoly_init(hsave, ctx);
        fmpz_mod_mpoly_init(gsave, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2 + 1, exp_bits2, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

            fmpz_mod_mpoly_mul(h, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            fmpz_mod_mpoly_set(hsave, h, ctx);
            fmpz_mod_mpoly_set(gsave, g, ctx);
            ret = fmpz_mod_mpoly_divides(k, h, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && fmpz_mod_mpoly_equal(k, f, ctx);

            if (!result)
            {
                flint_printf("FAIL: Check f*g/g = f sparse\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            if (   !fmpz_mod_mpoly_equal(h, hsave, ctx)
                || !fmpz_mod_mpoly_equal(g, gsave, ctx))
            {
                flint_printf("FAIL: Check f*g/g = f sparse input modification\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_set_num_threads(n_randint(state, max_threads) + 1);

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_clear(hsave, ctx);
        fmpz_mod_mpoly_clear(gsave, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f dense */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 6, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 150/n/n;
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
            fmpz_mod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bounds(g, state, len2 + 1, exp_bound2, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bounds(h, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mod_mpoly_mul(h, f, g, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            ret = fmpz_mod_mpoly_divides(k, h, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && fmpz_mod_mpoly_equal(k, f, ctx);

            if (!result)
            {
                flint_printf("FAIL: Check f*g/g = f dense\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check divisibility of random polys */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        fmpz * shifts, * strides;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 6, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        max_bound = 1 + 20/n;
        exp_bound = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound1 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        exp_bound2 = (mp_limb_t *) flint_malloc(ctx->minfo->nvars*sizeof(mp_limb_t));
        shifts = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        strides = (fmpz *) flint_malloc(ctx->minfo->nvars*sizeof(fmpz));
        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            exp_bound[j] = UWORD(1) << (FLINT_BITS - 1);
            exp_bound1[j] = n_randint(state, max_bound) + 1;
            exp_bound2[j] = n_randint(state, max_bound) + 1;
            fmpz_init(shifts + j);
            fmpz_init(strides + j);
            fmpz_randtest_unsigned(shifts + j, state, 100);
            fmpz_randtest_unsigned(strides + j, state, 100);
            fmpz_add_ui(strides + j, strides + j, 1);
        }

        for (j = 0; j < 4; j++)
        {
            fmpz_mod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bounds(g, state, len2 + 1, exp_bound2, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bounds(h, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpz_mod_mpoly_inflate(g, g, shifts, strides, ctx);

            ret = fmpz_mod_mpoly_divides(h, f, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            if (ret == 0)
                continue;
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            fmpz_mod_mpoly_mul(k, h, g, ctx);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            result = fmpz_mod_mpoly_equal(k, f, ctx);
            if (!result)
            {
                flint_printf("FAIL: Check divisibility of random polys\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            fmpz_clear(shifts + j);
            fmpz_clear(strides + j);
        }
        flint_free(shifts);
        flint_free(strides);

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f aliasing first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 6, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);

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
            fmpz_mod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bounds(g, state, len2 + 1, exp_bound2, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bounds(h, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mod_mpoly_set(h, f, ctx);
            fmpz_mod_mpoly_mul(h, h, g, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            ret = fmpz_mod_mpoly_divides(k, h, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && fmpz_mod_mpoly_equal(k, f, ctx);

            if (!result)
            {
                flint_printf("FAIL: Check f*g/g = f aliasing first\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f aliasing second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        slong n;

        fmpz_mod_mpoly_ctx_init_rand_bits_prime(ctx, state, 6, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(k, ctx);

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
            fmpz_mod_mpoly_randtest_bounds(f, state, len1, exp_bound1, ctx);
            fmpz_mod_mpoly_randtest_bounds(g, state, len2 + 1, exp_bound2, ctx);
            if (fmpz_mod_mpoly_is_zero(g, ctx))
                fmpz_mod_mpoly_one(g, ctx);
            fmpz_mod_mpoly_randtest_bounds(h, state, len, exp_bound, ctx);
            fmpz_mod_mpoly_randtest_bounds(k, state, len, exp_bound, ctx);

            fmpz_mod_mpoly_set(h, g, ctx);
            fmpz_mod_mpoly_mul(h, f, h, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            ret = fmpz_mod_mpoly_divides(k, h, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            fmpz_mod_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && fmpz_mod_mpoly_equal(k, f, ctx);

            if (!result)
            {
                flint_printf("FAIL: Check f*g/g = f aliasing second\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(k, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
