/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_quasidivrem_heap, state)
{
    int i, j, result;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, r1, q1;
        fmpz_t s1;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(q1, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_init(s1);

        len = n_randint(state, 100);
        len1 = n_randint(state, 10) + 50;
        len2 = n_randint(state, 10) + 50;

        exp_bits =  n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(q1, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(r1, state, len, coeff_bits, exp_bits, ctx);
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
            } while (g->length == 0);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_quasidivrem_heap(s1, q1, r1, h, g, ctx);
            fmpz_mpoly_assert_canonical(q1, ctx);
            fmpz_mpoly_assert_canonical(r1, ctx);
            fmpz_mpoly_remainder_strongtest(r1, g, ctx);

            result = fmpz_is_one(s1) && fmpz_mpoly_equal(q1, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni=%wd, j=%wd\n",i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(q1, ctx);
        fmpz_mpoly_clear(r1, ctx);
        fmpz_clear(s1);
    }

    /* Check aliasing of quotient with first argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, r1, t1, t2;
        fmpz_t s1, s2;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(s1);
        fmpz_init(s2);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound =  n_randint(state, 1 + 200/n/n) + 1;
        exp_bound1 = n_randint(state, 1 + 200/n/n) + 1;
        exp_bound2 = n_randint(state, 1 + 200/n/n) + 1;

        coeff_bits = n_randint(state, 50);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(r1, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);

            fmpz_mpoly_quasidivrem_heap(s1, h, r1, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_assert_canonical(r1, ctx);
            fmpz_mpoly_remainder_strongtest(r1, g, ctx);

            fmpz_mpoly_mul_johnson(t1, h, g, ctx);
            fmpz_mpoly_add(t1, t1, r1, ctx);
            fmpz_mpoly_scalar_mul_fmpz(t2, f, s1, ctx);

            fmpz_mpoly_quasidivrem_heap(s2, f, k, f, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            result = fmpz_mpoly_equal(t1, t2, ctx)
                    && fmpz_mpoly_equal(h, f, ctx)
                    && fmpz_mpoly_equal(r1, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                printf("Check aliasing of quotient with first argument\n");
                flint_printf("i = %wd  j = %wd\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(s1);
        fmpz_clear(s2);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(r1, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
    }

    /* Check aliasing of quotient with second argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, r1, t1, t2;
        fmpz_t s1, s2;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(s1);
        fmpz_init(s2);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound =  n_randint(state, 1 + 200/n/n) + 1;
        exp_bound1 = n_randint(state, 1 + 200/n/n) + 1;
        exp_bound2 = n_randint(state, 1 + 200/n/n) + 1;

        coeff_bits = n_randint(state, 50);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(r1, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);

            fmpz_mpoly_quasidivrem_heap(s1, h, r1, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_assert_canonical(r1, ctx);
            fmpz_mpoly_remainder_strongtest(r1, g, ctx);

            fmpz_mpoly_mul_johnson(t1, h, g, ctx);
            fmpz_mpoly_add(t1, t1, r1, ctx);
            fmpz_mpoly_scalar_mul_fmpz(t2, f, s1, ctx);

            fmpz_mpoly_quasidivrem_heap(s2, g, k, f, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            result = fmpz_mpoly_equal(t1, t2, ctx)
                    && fmpz_mpoly_equal(h, g, ctx)
                    && fmpz_mpoly_equal(r1, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                printf("Check aliasing of quotient with second argument\n");
                flint_printf("i = %wd  j = %wd\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(s1);
        fmpz_clear(s2);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(r1, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
    }

    /* Check aliasing of remainder with first argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, r1, t1, t2;
        fmpz_t s1, s2;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(s1);
        fmpz_init(s2);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound =  n_randint(state, 1 + 200/n/n) + 1;
        exp_bound1 = n_randint(state, 1 + 200/n/n) + 1;
        exp_bound2 = n_randint(state, 1 + 200/n/n) + 1;

        coeff_bits = n_randint(state, 50);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(r1, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);

            fmpz_mpoly_quasidivrem_heap(s1, h, r1, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_assert_canonical(r1, ctx);
            fmpz_mpoly_remainder_strongtest(r1, g, ctx);

            fmpz_mpoly_mul_johnson(t1, h, g, ctx);
            fmpz_mpoly_add(t1, t1, r1, ctx);
            fmpz_mpoly_scalar_mul_fmpz(t2, f, s1, ctx);

            fmpz_mpoly_quasidivrem_heap(s2, k, f, f, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);

            result = fmpz_mpoly_equal(t1, t2, ctx)
                    && fmpz_mpoly_equal(h, k, ctx)
                    && fmpz_mpoly_equal(r1, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                printf("Check aliasing of remainder with first argument\n");
                flint_printf("i = %wd  j = %wd\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(s1);
        fmpz_clear(s2);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(r1, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
    }

    /* Check aliasing of remainder with second argument */
    for (i = 0; i < 5 * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, k, r1, t1, t2;
        fmpz_t s1, s2;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpz_mpoly_ctx_init_rand(ctx, state, 10);

        fmpz_init(s1);
        fmpz_init(s2);

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(k, ctx);
        fmpz_mpoly_init(r1, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 10) + 1;

        n = FLINT_MAX(WORD(1), ctx->minfo->nvars);
        exp_bound =  n_randint(state, 1 + 200/n/n) + 1;
        exp_bound1 = n_randint(state, 1 + 200/n/n) + 1;
        exp_bound2 = n_randint(state, 1 + 200/n/n) + 1;

        coeff_bits = n_randint(state, 50);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpz_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (g->length == 0);
            fmpz_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);
            fmpz_mpoly_randtest_bound(r1, state, len, coeff_bits, exp_bound, ctx);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);

            fmpz_mpoly_quasidivrem_heap(s1, h, r1, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_assert_canonical(r1, ctx);
            fmpz_mpoly_remainder_strongtest(r1, g, ctx);

            fmpz_mpoly_mul_johnson(t1, h, g, ctx);
            fmpz_mpoly_add(t1, t1, r1, ctx);
            fmpz_mpoly_scalar_mul_fmpz(t2, f, s1, ctx);

            fmpz_mpoly_quasidivrem_heap(s2, k, g, f, g, ctx);
            fmpz_mpoly_assert_canonical(k, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);

            result = fmpz_mpoly_equal(t1, t2, ctx)
                    && fmpz_mpoly_equal(h, k, ctx)
                    && fmpz_mpoly_equal(r1, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                printf("Check aliasing of remainder with second argument\n");
                flint_printf("i = %wd  j = %wd\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(s1);
        fmpz_clear(s2);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(k, ctx);
        fmpz_mpoly_clear(r1, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
    }

    TEST_FUNCTION_END(state);
}
