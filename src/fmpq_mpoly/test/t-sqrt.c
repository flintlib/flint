/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpq_mpoly.h"

TEST_FUNCTION_START(fmpq_mpoly_sqrt, state)
{
    slong i, j, tmul = 10;

    /* Check sqrt(f^2) = +-f */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        flint_bitcnt_t coeff_bits, coeff_bits1;
        int sqr;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);

        exp_bits =  n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 100);
        coeff_bits1 = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits1, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpq_mpoly_mul(g, f, f, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);

            sqr = fmpq_mpoly_sqrt(h, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            if (!sqr)
            {
                flint_printf("FAIL: Check sqrt(f^2) returns 1\n");
                fflush(stdout);
                flint_abort();
            }

            if (!fmpq_mpoly_equal(h, f, ctx) &&
                !(fmpq_mpoly_neg(h, h, ctx), fmpq_mpoly_equal(h, f, ctx)))
            {
                flint_printf("FAIL: Check sqrt(f^2) = +-f\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check sqrt(random) */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        flint_bitcnt_t coeff_bits, coeff_bits1;
        int sqr;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100) + 1;

        exp_bits =  n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 100);
        coeff_bits1 = n_randint(state, 100) + 1;

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits1, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);

            sqr = fmpq_mpoly_sqrt(g, f, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);

            if (sqr)
            {
                fmpq_mpoly_mul(g, g, g, ctx);
                if (!fmpq_mpoly_equal(g, f, ctx))
                {
                    flint_printf("FAIL: Check sqrt(random)\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else if (!fmpq_mpoly_is_zero(g, ctx))
            {
               flint_printf("FAIL: Check nonsquare returns 0 sqrt\n");
               fflush(stdout);
               flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of square root with input */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        flint_bitcnt_t coeff_bits, coeff_bits1;
        int sqr1, sqr2;

        fmpq_mpoly_ctx_init_rand(ctx, state, 10);

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);

        exp_bits =  n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 100);
        coeff_bits1 = n_randint(state, 100);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits1, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpq_mpoly_mul(g, f, f, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);

            sqr1 = fmpq_mpoly_sqrt(h, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            sqr2 = fmpq_mpoly_sqrt(g, g, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);

            if (sqr1 != sqr2 || !fmpq_mpoly_equal(g, h, ctx))
            {
                printf("FAIL: Check aliasing\n");
                flint_printf("\n");
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
