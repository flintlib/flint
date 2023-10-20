/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_sqrt, state)
{
    slong i, j, tmul = 10;

    /* Check sqrt(f^2) = +-f */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        int sqr;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100) + 1;

        exp_bits =  n_randint(state, 5) + 1;
        exp_bits1 = n_randint(state, 5) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            fq_nmod_mpoly_mul(g, f, f, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);

            sqr = fq_nmod_mpoly_sqrt(h, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            if (!sqr)
            {
                flint_printf("FAIL: Check sqrt(f^2) returns 1\n");
                fflush(stdout);
                flint_abort();
            }

            if (!fq_nmod_mpoly_equal(h, f, ctx) &&
                !(fq_nmod_mpoly_neg(h, h, ctx), fq_nmod_mpoly_equal(h, f, ctx)))
            {
                flint_printf("FAIL: Check sqrt(f^2) = +-f\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check sqrt(random) */
    for (i = 0; i < 2*tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        int sqr;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100) + 1;

        exp_bits =  n_randint(state, 200) + 1;
        /* low bits: sqrt(random) is less reliable in positive char */
        exp_bits1 = n_randint(state, 20) + 1;

        for (j = 0; j < 10; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);

            sqr = fq_nmod_mpoly_sqrt(g, f, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);

            if (sqr)
            {
                fq_nmod_mpoly_mul(g, g, g, ctx);
                if (!fq_nmod_mpoly_equal(g, f, ctx))
                {
                    flint_printf("FAIL: Check sqrt(random)\n");
                    flint_printf("i = %wd, j = %wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
            else if (!fq_nmod_mpoly_is_zero(g, ctx))
            {
               flint_printf("FAIL: Check nonsquare returns 0 sqrt\n");
               fflush(stdout);
               flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing of square root with input */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        slong len, len1;
        flint_bitcnt_t exp_bits, exp_bits1;
        int sqr1, sqr2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 5);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);

        exp_bits =  n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            fq_nmod_mpoly_mul(g, f, f, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);

            sqr1 = fq_nmod_mpoly_sqrt(h, g, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            sqr2 = fq_nmod_mpoly_sqrt(g, g, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);

            if (sqr1 != sqr2 || !fq_nmod_mpoly_equal(g, h, ctx))
            {
                printf("FAIL: Check aliasing\n");
                flint_printf("\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
