/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_divrem, state)
{
    int i, j, result;

    /* Check f*g/g = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k, r;
        slong len, len1, len2;
        slong coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 50);
        len1 = n_randint(state, 50);
        len2 = n_randint(state, 50) + 1;

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            do {
                fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits + 1, exp_bits2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(k, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(r, state, len, coeff_bits, exp_bits, ctx);

            fmpq_mpoly_mul(h, f, g, ctx);

            fmpq_mpoly_divrem(k, r, h, g, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);
            fmpq_mpoly_assert_canonical(r, ctx);

            result = fmpq_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*g/g = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(k, ctx);
        fmpq_mpoly_clear(r, ctx);
    }

    /* Check f = g*q + r for random polys */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k, r;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        n = FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx));
        exp_bound = n_randint(state, 25/n + 1) + 2;
        exp_bound1 = n_randint(state, 35/n + 1) + 4;
        exp_bound2 = n_randint(state, 30/n + 1) + 2;

        coeff_bits = n_randint(state, 70);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));
            fmpq_mpoly_randtest_bound(h, state, len, coeff_bits, exp_bound, ctx);
            fmpq_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_divrem(h, r, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_assert_canonical(r, ctx);
            fmpq_mpoly_remainder_test(r, g, ctx);

            fmpq_mpoly_mul(k, h, g, ctx);
            fmpq_mpoly_add(k, k, r, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);

            result = fmpq_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
       }

       fmpq_mpoly_clear(f, ctx);
       fmpq_mpoly_clear(g, ctx);
       fmpq_mpoly_clear(h, ctx);
       fmpq_mpoly_clear(k, ctx);
       fmpq_mpoly_clear(r, ctx);
    }

    /* Check f = g*q + r for random polys, alias quotient and numerator */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k, r;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        n = FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx));
        exp_bound = n_randint(state, 35/n + 1) + 2;
        exp_bound1 = n_randint(state, 35/n + 1) + 4;
        exp_bound2 = n_randint(state, 30/n + 1) + 2;

        coeff_bits = n_randint(state, 70);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));

            fmpq_mpoly_set(h, f, ctx);
            fmpq_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_divrem(h, r, h, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_assert_canonical(r, ctx);
            fmpq_mpoly_remainder_test(r, g, ctx);

            fmpq_mpoly_mul(k, h, g, ctx);
            fmpq_mpoly_add(k, k, r, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);

            result = fmpq_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys, alias quotient and numerator\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
       }

       fmpq_mpoly_clear(f, ctx);
       fmpq_mpoly_clear(g, ctx);
       fmpq_mpoly_clear(h, ctx);
       fmpq_mpoly_clear(k, ctx);
       fmpq_mpoly_clear(r, ctx);
    }

    /* Check f = g*q + r for random polys, alias quotient and denominator */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k, r;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        n = FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx));
        exp_bound = n_randint(state, 25/n + 1) + 2;
        exp_bound1 = n_randint(state, 35/n + 1) + 4;
        exp_bound2 = n_randint(state, 30/n + 1) + 2;

        coeff_bits = n_randint(state, 70);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));

            fmpq_mpoly_set(h, g, ctx);
            fmpq_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_divrem(h, r, f, h, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_assert_canonical(r, ctx);
            fmpq_mpoly_remainder_test(r, g, ctx);

            fmpq_mpoly_mul(k, h, g, ctx);
            fmpq_mpoly_add(k, k, r, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);

            result = fmpq_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys, alias quotient and numerator\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
       }

       fmpq_mpoly_clear(f, ctx);
       fmpq_mpoly_clear(g, ctx);
       fmpq_mpoly_clear(h, ctx);
       fmpq_mpoly_clear(k, ctx);
       fmpq_mpoly_clear(r, ctx);
    }

    /* Check f = g*q + r for random polys, alias remainder and denominator */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k, r;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        n = FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx));
        exp_bound = n_randint(state, 25/n + 1) + 2;
        exp_bound1 = n_randint(state, 35/n + 1) + 4;
        exp_bound2 = n_randint(state, 30/n + 1) + 2;

        coeff_bits = n_randint(state, 70);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));

            fmpq_mpoly_set(r, g, ctx);
            fmpq_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_divrem(h, r, f, r, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_assert_canonical(r, ctx);
            fmpq_mpoly_remainder_test(r, g, ctx);

            fmpq_mpoly_mul(k, h, g, ctx);
            fmpq_mpoly_add(k, k, r, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);

            result = fmpq_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys, alias quotient and numerator\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
       }

       fmpq_mpoly_clear(f, ctx);
       fmpq_mpoly_clear(g, ctx);
       fmpq_mpoly_clear(h, ctx);
       fmpq_mpoly_clear(k, ctx);
       fmpq_mpoly_clear(r, ctx);
    }

    /* Check f = g*q + r for random polys, alias remainder and numerator */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k, r;
        slong len, len1, len2, exp_bound, exp_bound1, exp_bound2;
        slong coeff_bits;
        slong n;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k, ctx);
        fmpq_mpoly_init(r, ctx);

        len = n_randint(state, 20);
        len1 = n_randint(state, 20);
        len2 = n_randint(state, 20) + 1;

        n = FLINT_MAX(WORD(1), fmpq_mpoly_ctx_nvars(ctx));
        exp_bound = n_randint(state, 25/n + 1) + 2;
        exp_bound1 = n_randint(state, 35/n + 1) + 4;
        exp_bound2 = n_randint(state, 30/n + 1) + 2;

        coeff_bits = n_randint(state, 70);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bound(f, state, len1, coeff_bits, exp_bound1, ctx);
            do {
                fmpq_mpoly_randtest_bound(g, state, len2, coeff_bits + 1, exp_bound2, ctx);
            } while (fmpq_mpoly_is_zero(g, ctx));

            fmpq_mpoly_set(r, f, ctx);
            fmpq_mpoly_randtest_bound(k, state, len, coeff_bits, exp_bound, ctx);

            fmpq_mpoly_divrem(h, r, r, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_assert_canonical(r, ctx);
            fmpq_mpoly_remainder_test(r, g, ctx);

            fmpq_mpoly_mul(k, h, g, ctx);
            fmpq_mpoly_add(k, k, r, ctx);
            fmpq_mpoly_assert_canonical(k, ctx);

            result = fmpq_mpoly_equal(f, k, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f = g*q + r for random polys, alias quotient and numerator\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
       }

       fmpq_mpoly_clear(f, ctx);
       fmpq_mpoly_clear(g, ctx);
       fmpq_mpoly_clear(h, ctx);
       fmpq_mpoly_clear(k, ctx);
       fmpq_mpoly_clear(r, ctx);
    }

    TEST_FUNCTION_END(state);
}
