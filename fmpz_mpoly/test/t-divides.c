/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mpoly.h"
#include "nmod_mpoly.h"

int
main(void)
{
    int i, j, result, ret, max_threads = 5, tmul = 25;
    FLINT_TEST_INIT(state);

    flint_printf("divides....");
    fflush(stdout);

    /* Check f*g/g = f sparse */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, hsave, gsave;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2, coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 6);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(hsave, ctx);
        fmpz_mpoly_init(gsave, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_mul(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_set(hsave, h, ctx);
            fmpz_mpoly_set(gsave, g, ctx);
            ret = fmpz_mpoly_divides(k, h, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            fmpz_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && fmpz_mpoly_equal(k, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f sparse\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                flint_abort();
            }

            if (   !fmpz_mpoly_equal(h, hsave, ctx)
                || !fmpz_mpoly_equal(g, gsave, ctx))
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f sparse input modification\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(hsave, ctx);
        fmpz_mpoly_clear(gsave, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check divisibility of random polys */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2, coeff_bits;
        fmpz * shifts, * strides;

        fmpz_mpoly_ctx_init_rand(ctx, state, 6);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50);

        coeff_bits = n_randint(state, 200) + 1;

        max_bound = 1 + 20/ctx->minfo->nvars;
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
            fmpz_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bounds(g, state, len2 + 1, coeff_bits, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bounds(k, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_inflate(f, f, shifts, strides, ctx);
            fmpz_mpoly_inflate(g, g, shifts, strides, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            ret = fmpz_mpoly_divides(h, f, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            if (ret == 0)
                continue;
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_mul(k, h, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            result = fmpz_mpoly_equal(k, f, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check divisibility of random polys\n"
                                                   "i = %wd, j = %wd\n", i ,j);
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

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f aliasing first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        flint_bitcnt_t coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 6);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200) + 1;

        max_bound = 1 + 100/ctx->minfo->nvars/ctx->minfo->nvars;
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
            fmpz_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bounds(g, state, len2 + 1, coeff_bits, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bounds(k, state, len, coeff_bits, exp_bound, ctx);

            flint_set_num_threads(n_randint(state, max_threads) + 1);

            fmpz_mpoly_set(h, f, ctx);
            fmpz_mpoly_mul(h, h, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            ret = fmpz_mpoly_divides(k, h, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            fmpz_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && fmpz_mpoly_equal(k, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing first\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    /* Check f*g/g = f aliasing second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k;
        slong len, len1, len2;
        mp_limb_t max_bound, * exp_bound, * exp_bound1, * exp_bound2;
        flint_bitcnt_t coeff_bits;

        fmpz_mpoly_ctx_init_rand(ctx, state, 6);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        coeff_bits = n_randint(state, 200) + 1;

        max_bound = 1 + 100/ctx->minfo->nvars/ctx->minfo->nvars;
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
            fmpz_mpoly_randtest_bounds(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bounds(g, state, len2 + 1, coeff_bits, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bounds(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bounds(k, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_set(h, g, ctx);
            fmpz_mpoly_mul(h, f, h, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            ret = fmpz_mpoly_divides(k, h, g, ctx);
            FLINT_ASSERT(ret == 0 || ret == 1);
            fmpz_mpoly_assert_canonical(k, ctx);
            result = (ret == 1) && fmpz_mpoly_equal(k, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f aliasing second\n"
                                                   "i = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        flint_free(exp_bound);
        flint_free(exp_bound1);
        flint_free(exp_bound2);

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

