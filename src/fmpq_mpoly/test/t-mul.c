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

TEST_FUNCTION_START(fmpq_mpoly_mul, state)
{
    int i, j, result;

    /* Check f*(g + h) = f*g + f*h */
    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, k1, k2, t1, t2;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(k1, ctx);
        fmpq_mpoly_init(k2, ctx);
        fmpq_mpoly_init(t1, ctx);
        fmpq_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        fmpq_mpoly_randtest_bits(k1, state, len, coeff_bits, exp_bits, ctx);
        fmpq_mpoly_randtest_bits(k2, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_randtest_bits(h, state, len2, coeff_bits, exp_bits2, ctx);

            fmpq_mpoly_add(t1, g, h, ctx);
            fmpq_mpoly_assert_canonical(t1, ctx);
            fmpq_mpoly_mul(k1, f, t1, ctx);
            fmpq_mpoly_assert_canonical(k1, ctx);
            fmpq_mpoly_mul(t1, f, g, ctx);
            fmpq_mpoly_assert_canonical(t1, ctx);
            fmpq_mpoly_mul(t2, f, h, ctx);
            fmpq_mpoly_assert_canonical(t2, ctx);
            fmpq_mpoly_add(k2, t1, t2, ctx);
            fmpq_mpoly_assert_canonical(k2, ctx);
            result = fmpq_mpoly_equal(k1, k2, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check f*(g + h) = f*g + f*h\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(k1, ctx);
        fmpq_mpoly_clear(k2, ctx);
        fmpq_mpoly_clear(t1, ctx);
        fmpq_mpoly_clear(t2, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing first argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpq_mpoly_mul(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_mul(f, f, g, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);
            result = fmpq_mpoly_equal(h, f, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing second argument */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        slong len, len1, len2;
        flint_bitcnt_t coeff_bits, exp_bits, exp_bits1, exp_bits2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 20);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpq_mpoly_mul(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_mul(g, f, g, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);
            result = fmpq_mpoly_equal(h, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing second arg\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
