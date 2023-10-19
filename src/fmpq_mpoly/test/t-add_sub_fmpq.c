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

TEST_FUNCTION_START(fmpq_mpoly_add_sub_fmpq, state)
{
    slong i, j;

    /* Check (f + a) - a = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h;
        fmpq_t c;
        slong len, coeff_bits, exp_bits;

        fmpq_init(c);
        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 10; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpq_randtest(c, state, n_randint(state, 200) + 1);

            fmpq_mpoly_add_fmpq(g, f, c, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);

            fmpq_mpoly_sub_fmpq(h, g, c, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            if (!fmpq_mpoly_equal(f, h, ctx))
            {
                flint_printf("FAIL: Check (f + a) - a = f\n");
                fflush(stdout);
                flint_abort();
            }

            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(g, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);

            fmpq_randtest(c, state, n_randint(state, 200) + 1);

            fmpq_mpoly_sub_fmpq(g, f, c, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);

            fmpq_mpoly_add_fmpq(h, g, c, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            if (!fmpq_mpoly_equal(f, h, ctx))
            {
                flint_printf("FAIL: Check (f - a) + a = f\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);

        fmpq_clear(c);
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g;
        fmpq_t c;
        slong len, coeff_bits, exp_bits;

        fmpq_init(c);
        fmpq_mpoly_ctx_init_rand(ctx, state, 20);
        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;
        coeff_bits = n_randint(state, 200);

        for (j = 0; j < 10; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
            fmpq_mpoly_set(g, f, ctx);

            fmpq_randtest(c, state, n_randint(state, 200) + 1);

            fmpq_mpoly_add_fmpq(f, f, c, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);

            fmpq_mpoly_sub_fmpq(f, f, c, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);

            if (!fmpq_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL: check aliasing\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);

        fmpq_clear(c);
    }

    TEST_FUNCTION_END(state);
}
