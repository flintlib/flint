/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fq_nmod_mpoly.h"

int
main(void)
{
    int result;
    slong i, j, tmul = 10;
    FLINT_TEST_INIT(state);

    flint_printf("divrem_monagan_pearce....");
    fflush(stdout);

    /* Check f*g/g = f */
    for (i = 0; i < 5 * tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k, r;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);
        fq_nmod_mpoly_init(r, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bits2, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(r, state, len, exp_bits, ctx);

            fq_nmod_mpoly_mul_johnson(h, f, g, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(k, r, h, g, ctx);
            fq_nmod_mpoly_assert_canonical(k, ctx);
            fq_nmod_mpoly_assert_canonical(r, ctx);

            result = fq_nmod_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);  
        fq_nmod_mpoly_clear(g, ctx);  
        fq_nmod_mpoly_clear(h, ctx);  
        fq_nmod_mpoly_clear(k, ctx);  
        fq_nmod_mpoly_clear(r, ctx);  

        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check f = g*q + r for random polys */
    for (i = 0; i < 5 * tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k, r;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);
        fq_nmod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 15);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound1 = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound2 = n_randint(state, 50/ctx->minfo->nvars) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(h, r, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_assert_canonical(r, ctx);
            fq_nmod_mpoly_remainder_strongtest(r, g, ctx);

            fq_nmod_mpoly_mul_johnson(k, h, g, ctx);
            fq_nmod_mpoly_add(k, k, r, ctx);
            fq_nmod_mpoly_assert_canonical(k, ctx);

            result = fq_nmod_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);  
        fq_nmod_mpoly_clear(g, ctx);  
        fq_nmod_mpoly_clear(h, ctx);  
        fq_nmod_mpoly_clear(k, ctx);  
        fq_nmod_mpoly_clear(r, ctx);  
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of quotient with first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, r1, r2;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(r1, ctx);
        fq_nmod_mpoly_init(r2, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound1 = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound2 = n_randint(state, 50/ctx->minfo->nvars) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(r1, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(r2, state, len, exp_bound, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_assert_canonical(r1, ctx);
            fq_nmod_mpoly_remainder_strongtest(r1, g, ctx);
            fq_nmod_mpoly_divrem_monagan_pearce(f, r2, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_assert_canonical(r2, ctx);
            fq_nmod_mpoly_remainder_strongtest(r2, g, ctx);

            result =    fq_nmod_mpoly_equal(h, f, ctx)
                     && fq_nmod_mpoly_equal(r1, r2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of quotient with first argument\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);  
        fq_nmod_mpoly_clear(g, ctx);  
        fq_nmod_mpoly_clear(h, ctx);  
        fq_nmod_mpoly_clear(r1, ctx);  
        fq_nmod_mpoly_clear(r2, ctx);  
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of quotient with second argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, r1, r2;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(r1, ctx);
        fq_nmod_mpoly_init(r2, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound1 = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound2 = n_randint(state, 50/ctx->minfo->nvars) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(r1, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(r2, state, len, exp_bound, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_assert_canonical(r1, ctx);
            fq_nmod_mpoly_remainder_strongtest(r1, g, ctx);
            fq_nmod_mpoly_divrem_monagan_pearce(g, r2, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_assert_canonical(r2, ctx);

            result =    fq_nmod_mpoly_equal(h, g, ctx)
                     && fq_nmod_mpoly_equal(r1, r2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of quotient with second argument\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);  
        fq_nmod_mpoly_clear(g, ctx);  
        fq_nmod_mpoly_clear(h, ctx);  
        fq_nmod_mpoly_clear(r1, ctx);  
        fq_nmod_mpoly_clear(r2, ctx);  
        fq_nmod_mpoly_ctx_clear(ctx);
    }
	
    /* Check aliasing of remainder with first argument */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k, r1;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);
        fq_nmod_mpoly_init(r1, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 15);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound1 = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound2 = n_randint(state, 50/ctx->minfo->nvars) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(r1, state, len, exp_bound, ctx);

            fq_nmod_mpoly_mul_johnson(h, f, g, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_assert_canonical(r1, ctx);
            fq_nmod_mpoly_remainder_strongtest(r1, g, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(k, f, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(k, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);
            fq_nmod_mpoly_remainder_strongtest(f, g, ctx);

            result =    fq_nmod_mpoly_equal(h, k, ctx)
                     && fq_nmod_mpoly_equal(r1, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of remainder with first argument\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);  
        fq_nmod_mpoly_clear(g, ctx);  
        fq_nmod_mpoly_clear(h, ctx);  
        fq_nmod_mpoly_clear(k, ctx);  
        fq_nmod_mpoly_clear(r1, ctx);  
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of remainder with second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k, r1;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);
        fq_nmod_mpoly_init(r1, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 15);
        len2 = n_randint(state, 10) + 1;

        exp_bound = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound1 = n_randint(state, 50/ctx->minfo->nvars) + 1;
        exp_bound2 = n_randint(state, 50/ctx->minfo->nvars) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            fq_nmod_mpoly_randtest_bound(g, state, len2, exp_bound2, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bound(h, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(k, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(r1, state, len, exp_bound, ctx);

            fq_nmod_mpoly_mul_johnson(h, f, g, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(h, r1, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_assert_canonical(r1, ctx);
            fq_nmod_mpoly_remainder_strongtest(r1, g, ctx);

            fq_nmod_mpoly_divrem_monagan_pearce(k, g, f, g, ctx);
            fq_nmod_mpoly_assert_canonical(k, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);

            result =    fq_nmod_mpoly_equal(h, k, ctx)
                     && fq_nmod_mpoly_equal(r1, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing of remainder with second argument\ni = %wd, j = %wd\n", i ,j);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);  
        fq_nmod_mpoly_clear(g, ctx);  
        fq_nmod_mpoly_clear(h, ctx);  
        fq_nmod_mpoly_clear(k, ctx);  
        fq_nmod_mpoly_clear(r1, ctx);  
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
