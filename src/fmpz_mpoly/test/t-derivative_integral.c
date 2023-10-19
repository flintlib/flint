/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mpoly.h"

TEST_FUNCTION_START(fmpz_mpoly_derivative_integral, state)
{
    int i, j, result;
    slong tmul = 5;

    /* randomized testing doesn't catch exponent overflow in integral */
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g;
        fmpz_t s;
        const char* vars[] = {"x","y","z","w","s","t","u","v"};

        fmpz_mpoly_ctx_init(ctx, 7, ORD_DEGLEX);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_init(s);

        result = fmpz_mpoly_set_str_pretty(f, "-2*x^2*(y+z)+z*y^126", vars, ctx);
        if (result)
        {
            printf("FAIL\n");
            flint_printf("set_str_pretty\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_mpoly_integral(g, s, f, 1, ctx);
        fmpz_mpoly_derivative(g, g, 1, ctx);
        fmpz_mpoly_scalar_divexact_fmpz(g, g, s, ctx);

        result = fmpz_mpoly_equal(f, g, ctx);
        if (!result)
        {
            printf("FAIL\n");
            flint_printf("manual integral check\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(s);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
    }

    /* Check d(f*g) = df*g + f*dg */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, fp, gp, hp, t1, t2;
        slong nvars, len, len1, len2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;

        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);
        fmpz_mpoly_init(hp, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest_bits(hp, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_randtest_bits(fp, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_randtest_bits(gp, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            idx = n_randint(state, nvars);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_derivative(hp, h, idx, ctx);
            fmpz_mpoly_assert_canonical(hp, ctx);
            fmpz_mpoly_derivative(fp, f, idx, ctx);
            fmpz_mpoly_assert_canonical(fp, ctx);
            fmpz_mpoly_derivative(gp, g, idx, ctx);
            fmpz_mpoly_assert_canonical(gp, ctx);
            fmpz_mpoly_mul_johnson(t1, f, gp, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            fmpz_mpoly_mul_johnson(t2, g, fp, ctx);
            fmpz_mpoly_assert_canonical(t2, ctx);
            fmpz_mpoly_add(t1, t1, t2, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);

            result = fmpz_mpoly_equal(hp, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg\ni = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(h, ctx);
        fmpz_mpoly_clear(fp, ctx);
        fmpz_mpoly_clear(gp, ctx);
        fmpz_mpoly_clear(hp, ctx);
        fmpz_mpoly_clear(t1, ctx);
        fmpz_mpoly_clear(t2, ctx);
    }

    /* Check d(f*g) = df*g + f*dg with aliasing */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, h, fp, gp, t1, t2;
        slong nvars, len, len1, len2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;

        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(h, ctx);
        fmpz_mpoly_init(fp, ctx);
        fmpz_mpoly_init(gp, ctx);
        fmpz_mpoly_init(t1, ctx);
        fmpz_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 200);

        fmpz_mpoly_randtest_bits(h, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_randtest_bits(fp, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_randtest_bits(gp, state, len, coeff_bits, exp_bits, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);

            idx = n_randint(state, nvars);

            fmpz_mpoly_mul_johnson(h, f, g, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);

            fmpz_mpoly_derivative(h, h, idx, ctx);
            fmpz_mpoly_assert_canonical(h, ctx);
            fmpz_mpoly_set(fp, f, ctx);
            fmpz_mpoly_derivative(fp, fp, idx, ctx);
            fmpz_mpoly_assert_canonical(fp, ctx);
            fmpz_mpoly_set(gp, g, ctx);
            fmpz_mpoly_derivative(gp, gp, idx, ctx);
            fmpz_mpoly_assert_canonical(gp, ctx);

            fmpz_mpoly_mul_johnson(t1, f, gp, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);
            fmpz_mpoly_mul_johnson(t2, g, fp, ctx);
            fmpz_mpoly_assert_canonical(t2, ctx);
            fmpz_mpoly_add(t1, t1, t2, ctx);
            fmpz_mpoly_assert_canonical(t1, ctx);

            result = fmpz_mpoly_equal(h, t1, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(f*g) = df*g + f*dg with aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
           }
       }

       fmpz_mpoly_clear(f, ctx);
       fmpz_mpoly_clear(g, ctx);
       fmpz_mpoly_clear(h, ctx);
       fmpz_mpoly_clear(fp, ctx);
       fmpz_mpoly_clear(gp, ctx);
       fmpz_mpoly_clear(t1, ctx);
       fmpz_mpoly_clear(t2, ctx);
    }

    /* Check d(int(f)) = f with aliasing */
    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t f, g, f1;
        slong nvars, len, len1, len2;
        slong idx, coeff_bits, exp_bits, exp_bits1, exp_bits2;
        fmpz_t s, s1;

        nvars = n_randint(state, 20) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, mpoly_ordering_randtest(state));

        fmpz_init(s);
        fmpz_init(s1);
        fmpz_mpoly_init(f, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(f1, ctx);

        len = n_randint(state, 10);
        len1 = n_randint(state, 10);
        len2 = n_randint(state, 10);

        exp_bits = n_randint(state, 200) + 1;
        exp_bits1 = n_randint(state, 200) + 1;
        exp_bits2 = n_randint(state, 200) + 1;

        coeff_bits = n_randint(state, 20);

        fmpz_mpoly_randtest_bits(f, state, len, coeff_bits, exp_bits, ctx);
        fmpz_mpoly_randtest_bits(g, state, len2, coeff_bits, exp_bits2, ctx);
        fmpz_mpoly_randtest_bits(f1, state, len1, coeff_bits, exp_bits1, ctx);

        for (j = 0; j < 4; j++)
        {
            fmpz_mpoly_randtest_bits(f, state, len1, coeff_bits, exp_bits1, ctx);
            idx = n_randint(state, nvars);

            fmpz_mpoly_set(g, f, ctx);
            fmpz_mpoly_assert_canonical(g, ctx);
            fmpz_mpoly_integral(f1, s1, f, idx, ctx);
            fmpz_mpoly_assert_canonical(f1, ctx);
            fmpz_mpoly_integral(f, s, f, idx, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            result = fmpz_equal(s, s1) && fmpz_mpoly_equal(f, f1, ctx);
            fmpz_mpoly_derivative(f, f, idx, ctx);
            fmpz_mpoly_assert_canonical(f, ctx);
            fmpz_mpoly_scalar_mul_fmpz(g, g, s, ctx);
            result = result && fmpz_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check d(int(f)) = f with aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(s);
        fmpz_clear(s1);
        fmpz_mpoly_clear(f, ctx);
        fmpz_mpoly_clear(g, ctx);
        fmpz_mpoly_clear(f1, ctx);
    }

    TEST_FUNCTION_END(state);
}
