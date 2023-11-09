/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_divrem_ideal_monagan_pearce, state)
{
    int result;
    slong i, j, w;

    /* Check f*g/g = f */
    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k, r;
        slong len, len1, len2;
        flint_bitcnt_t exp_bits, exp_bits1, exp_bits2;
        fq_nmod_mpoly_struct * qarr[1], * darr[1];

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);
        fq_nmod_mpoly_init(r, ctx);

        len = n_randint(state, 30);
        len1 = n_randint(state, 30);
        len2 = n_randint(state, 30) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2 + 1, ctx);
            if (fq_nmod_mpoly_is_zero(g, ctx))
                fq_nmod_mpoly_one(g, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(r, state, len, exp_bits, ctx);

            fq_nmod_mpoly_mul_johnson(h, f, g, ctx);

            qarr[0] = k;
            darr[0] = g;

            fq_nmod_mpoly_divrem_ideal_monagan_pearce(qarr, r, h, darr, 1, ctx);
            fq_nmod_mpoly_assert_canonical(qarr[0], ctx);
            fq_nmod_mpoly_assert_canonical(r, ctx);

            result = fq_nmod_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
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

    /* Check f = g1*q1 + ... + gn*qn + r for random polys */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, r, k1, k2;
        fq_nmod_mpoly_struct * g, * q;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        fq_nmod_mpoly_struct * qarr[5], * darr[5];
        slong n;

        num = n_randint(state, 5) + 1;
        g = (fq_nmod_mpoly_struct *) flint_malloc(num*sizeof(fq_nmod_mpoly_struct));
        q = (fq_nmod_mpoly_struct *) flint_malloc(num*sizeof(fq_nmod_mpoly_struct));

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        for (w = 0; w < num; w++)
        {
            fq_nmod_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            fq_nmod_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(k1, ctx);
        fq_nmod_mpoly_init(k2, ctx);
        fq_nmod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 10/n + 1) + 2;
        exp_bound1 = n_randint(state, 25/n + 1) + 2;
        exp_bound2 = n_randint(state, 20/n + 1) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            for (w = 0; w < num; w++)
            {
                fq_nmod_mpoly_randtest_bound(darr[w], state, len2, exp_bound2 + 1, ctx);
                if (fq_nmod_mpoly_is_zero(darr[w], ctx))
                    fq_nmod_mpoly_one(darr[w], ctx);
                fq_nmod_mpoly_randtest_bound(qarr[w], state, len, exp_bound, ctx);
            }
            fq_nmod_mpoly_randtest_bound(k1, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(k2, state, len, exp_bound, ctx);

            fq_nmod_mpoly_divrem_ideal_monagan_pearce(qarr, r, f, darr, num, ctx);
            fq_nmod_mpoly_assert_canonical(r, ctx);
            for (w = 0; w < num; w++)
            {
                fq_nmod_mpoly_assert_canonical(qarr[w], ctx);
                fq_nmod_mpoly_remainder_strongtest(r, darr[w], ctx);
            }

            fq_nmod_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                fq_nmod_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fq_nmod_mpoly_add(k2, k2, k1, ctx);
            }
            fq_nmod_mpoly_add(k2, k2, r, ctx);

            result = fq_nmod_mpoly_equal(f, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g1*q1 + ... + gn*qn + r for random polys\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            fq_nmod_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fq_nmod_mpoly_clear(darr[w], ctx);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(k1, ctx);
        fq_nmod_mpoly_clear(k2, ctx);
        fq_nmod_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, r, k1, k2;
        fq_nmod_mpoly_struct * g, * q;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2, num;
        fq_nmod_mpoly_struct * qarr[5], * darr[5];
        slong n;

        num = n_randint(state, 5) + 1;
        g = (fq_nmod_mpoly_struct *) flint_malloc(num*sizeof(fq_nmod_mpoly_struct));
        q = (fq_nmod_mpoly_struct *) flint_malloc(num*sizeof(fq_nmod_mpoly_struct));

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        for (w = 0; w < num; w++)
        {
            fq_nmod_mpoly_init(g + w, ctx);
            darr[w] = g + w;

            fq_nmod_mpoly_init(q + w, ctx);
            qarr[w] = q + w;
        }

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(k1, ctx);
        fq_nmod_mpoly_init(k2, ctx);
        fq_nmod_mpoly_init(r, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound = n_randint(state, 10/n + 1) + 2;
        exp_bound1 = n_randint(state, 25/n + 1) + 2;
        exp_bound2 = n_randint(state, 20/n + 1) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bound(f, state, len1, exp_bound1, ctx);
            for (w = 0; w < num; w++)
            {
                fq_nmod_mpoly_randtest_bound(darr[w], state, len2, exp_bound2 + 1, ctx);
                if (fq_nmod_mpoly_is_zero(darr[w], ctx))
                    fq_nmod_mpoly_one(darr[w], ctx);
                fq_nmod_mpoly_randtest_bound(qarr[w], state, len, exp_bound, ctx);
            }
            fq_nmod_mpoly_randtest_bound(k1, state, len, exp_bound, ctx);
            fq_nmod_mpoly_randtest_bound(k2, state, len, exp_bound, ctx);

            fq_nmod_mpoly_set(r, f, ctx);

            fq_nmod_mpoly_divrem_ideal_monagan_pearce(qarr, f, f, darr, num, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);
            for (w = 0; w < num; w++)
            {
                fq_nmod_mpoly_assert_canonical(qarr[w], ctx);
                fq_nmod_mpoly_remainder_strongtest(f, darr[w], ctx);
            }

            fq_nmod_mpoly_zero(k2, ctx);
            for (w = 0; w < num; w++)
            {
                fq_nmod_mpoly_mul_johnson(k1, qarr[w], darr[w], ctx);
                fq_nmod_mpoly_add(k2, k2, k1, ctx);
            }
            fq_nmod_mpoly_add(k2, k2, f, ctx);

            result = fq_nmod_mpoly_equal(r, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        for (w = 0; w < num; w++)
            fq_nmod_mpoly_clear(qarr[w], ctx);
        for (w = 0; w < num; w++)
            fq_nmod_mpoly_clear(darr[w], ctx);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(k1, ctx);
        fq_nmod_mpoly_clear(k2, ctx);
        fq_nmod_mpoly_clear(r, ctx);

        flint_free(g);
        flint_free(q);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
