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

TEST_FUNCTION_START(fmpq_mpoly_derivative_integral, state)
{
    int i, j, result;

    /* Check d(f*g) = df*g + f*dg */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, fp, gp, hp, t1, t2;
        slong nvars, len, len1, len2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;

        nvars = n_randint(state, 20) + 1;

        fmpq_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(fp, ctx);
        fmpq_mpoly_init(gp, ctx);
        fmpq_mpoly_init(hp, ctx);
        fmpq_mpoly_init(t1, ctx);
        fmpq_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);

        fmpq_mpoly_randtest_bits(hp, state, len, coeff_bits, exp_bits, ctx);
        fmpq_mpoly_randtest_bits(fp, state, len, coeff_bits, exp_bits, ctx);
        fmpq_mpoly_randtest_bits(gp, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            idx = n_randint(state, nvars);

            fmpq_mpoly_mul(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            fmpq_mpoly_derivative(hp, h, idx, ctx);
            fmpq_mpoly_assert_canonical(hp, ctx);
            fmpq_mpoly_derivative(fp, f, idx, ctx);
            fmpq_mpoly_assert_canonical(fp, ctx);
            fmpq_mpoly_derivative(gp, g, idx, ctx);
            fmpq_mpoly_assert_canonical(gp, ctx);

            fmpq_mpoly_mul(t1, f, gp, ctx);
            fmpq_mpoly_assert_canonical(t1, ctx);
            fmpq_mpoly_mul(t2, g, fp, ctx);
            fmpq_mpoly_assert_canonical(t2, ctx);
            fmpq_mpoly_add(t1, t1, t2, ctx);
            fmpq_mpoly_assert_canonical(t1, ctx);

            result = fmpq_mpoly_equal(hp, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(h, ctx);
        fmpq_mpoly_clear(fp, ctx);
        fmpq_mpoly_clear(gp, ctx);
        fmpq_mpoly_clear(hp, ctx);
        fmpq_mpoly_clear(t1, ctx);
        fmpq_mpoly_clear(t2, ctx);
    }

    /* Check d(f*g) = df*g + f*dg with aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, h, fp, gp, t1, t2;
        slong nvars, len, len1, len2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;

        nvars = n_randint(state, 20) + 1;

        fmpq_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(h, ctx);
        fmpq_mpoly_init(fp, ctx);
        fmpq_mpoly_init(gp, ctx);
        fmpq_mpoly_init(t1, ctx);
        fmpq_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);

        fmpq_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
        fmpq_mpoly_randtest_bits(fp, state, len, coeff_bits, exp_bits, ctx);
        fmpq_mpoly_randtest_bits(gp, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            idx = n_randint(state, nvars);

            fmpq_mpoly_mul(h, f, g, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);

            fmpq_mpoly_derivative(h, h, idx, ctx);
            fmpq_mpoly_assert_canonical(h, ctx);
            fmpq_mpoly_set(fp, f, ctx);
            fmpq_mpoly_derivative(fp, fp, idx, ctx);
            fmpq_mpoly_assert_canonical(fp, ctx);
            fmpq_mpoly_set(gp, g, ctx);
            fmpq_mpoly_derivative(gp, gp, idx, ctx);
            fmpq_mpoly_assert_canonical(gp, ctx);

            fmpq_mpoly_mul(t1, f, gp, ctx);
            fmpq_mpoly_assert_canonical(t1, ctx);
            fmpq_mpoly_mul(t2, g, fp, ctx);
            fmpq_mpoly_assert_canonical(t2, ctx);
            fmpq_mpoly_add(t1, t1, t2, ctx);
            fmpq_mpoly_assert_canonical(t1, ctx);

            result = fmpq_mpoly_equal(h, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg with aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
       }

       fmpq_mpoly_clear(f, ctx);
       fmpq_mpoly_clear(g, ctx);
       fmpq_mpoly_clear(h, ctx);
       fmpq_mpoly_clear(fp, ctx);
       fmpq_mpoly_clear(gp, ctx);
       fmpq_mpoly_clear(t1, ctx);
       fmpq_mpoly_clear(t2, ctx);
    }

    /* Check d(int(f)) = f with aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpq_mpoly_ctx_t ctx;
        fmpq_mpoly_t f, g, f1;
        slong nvars, len, len1, len2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;

        nvars = n_randint(state, 20) + 1;

        fmpq_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        fmpq_mpoly_init(f, ctx);
        fmpq_mpoly_init(g, ctx);
        fmpq_mpoly_init(f1, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);

        fmpq_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
        fmpq_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
        fmpq_mpoly_randtest_bits(f1, state, len1, coeff_bits, exp_bits1, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpq_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            idx = n_randint(state, nvars);

            fmpq_mpoly_set(g, f, ctx);
            fmpq_mpoly_assert_canonical(g, ctx);
            fmpq_mpoly_integral(f1, f, idx, ctx);
            fmpq_mpoly_assert_canonical(f1, ctx);
            fmpq_mpoly_integral(f, f, idx, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);
            result = fmpq_mpoly_equal(f, f1, ctx);
            fmpq_mpoly_derivative(f, f, idx, ctx);
            fmpq_mpoly_assert_canonical(f, ctx);
            result = result && fmpq_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(int(f)) = f with aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpq_mpoly_clear(f, ctx);
        fmpq_mpoly_clear(g, ctx);
        fmpq_mpoly_clear(f1, ctx);
    }

    TEST_FUNCTION_END(state);
}
