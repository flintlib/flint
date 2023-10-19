/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_scalar_addmul_fmpz, state)
{
    slong i, j;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h, t1, t2;
        fmpz_t a, b;
        slong len, exp_bits;

        fmpz_init(a);
        fmpz_init(b);

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);

        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_mod_mpoly_init(t1, ctx);
        fmpz_mod_mpoly_init(t2, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;

        for (j = 0; j < 10; j++)
        {
            fmpz_one(a);
            fmpz_randtest(b, state, n_randint(state, 200));

            fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(t1, state, len, exp_bits, ctx);

            len = n_randint(state, 100);
            exp_bits = n_randint(state, 200) + 1;

            fmpz_mod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fmpz_mod_mpoly_randtest_bits(t2, state, len, exp_bits, ctx);

            fmpz_mod_mpoly_scalar_addmul_fmpz(f, g, h, b, ctx);
            fmpz_mod_mpoly_assert_canonical(f, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(t1, g, a, ctx);
            fmpz_mod_mpoly_scalar_mul_fmpz(t2, h, b, ctx);
            fmpz_mod_mpoly_add(t1, t1, t2, ctx);
            if (!fmpz_mod_mpoly_equal(f, t1, ctx))
            {
                flint_printf("FAIL: check addmul definition\n");
                flint_printf("i = %wd, j = %wd\n", i,j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_set(t1, g, ctx);
            fmpz_mod_mpoly_scalar_addmul_fmpz(t1, t1, h, b, ctx);
            fmpz_mod_mpoly_assert_canonical(t1, ctx);
            if (!fmpz_mod_mpoly_equal(f, t1, ctx))
            {
                flint_printf("FAIL: check aliasing first argument\n");
                flint_printf("i = %wd, j = %wd\n", i,j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_set(t1, h, ctx);
            fmpz_mod_mpoly_scalar_addmul_fmpz(t1, g, t1, b, ctx);
            fmpz_mod_mpoly_assert_canonical(t1, ctx);
            if (!fmpz_mod_mpoly_equal(f, t1, ctx))
            {
                flint_printf("FAIL: check aliasing second argument\n");
                flint_printf("i = %wd, j = %wd\n", i,j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_scalar_addmul_fmpz(f, g, g, b, ctx);
            fmpz_mod_mpoly_set(t1, g, ctx);
            fmpz_mod_mpoly_scalar_addmul_fmpz(t1, t1, t1, b, ctx);
            fmpz_mod_mpoly_assert_canonical(t1, ctx);
            if (!fmpz_mod_mpoly_equal(f, t1, ctx))
            {
                flint_printf("FAIL: check aliasing both arguments\n");
                flint_printf("i = %wd, j = %wd\n", i,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_clear(t1, ctx);
        fmpz_mod_mpoly_clear(t2, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);

        fmpz_clear(a);
        fmpz_clear(b);
    }

    TEST_FUNCTION_END(state);
}
