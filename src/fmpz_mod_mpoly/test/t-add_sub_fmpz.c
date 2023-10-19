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

TEST_FUNCTION_START(fmpz_mod_mpoly_add_sub_fmpz, state)
{
    slong i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, g, h;
        fmpz_t c;
        slong len, len1, len2;
        slong exp_bits, exp_bits1, exp_bits2;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        fmpz_mod_mpoly_init(f, ctx);
        fmpz_mod_mpoly_init(g, ctx);
        fmpz_mod_mpoly_init(h, ctx);
        fmpz_init(c);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fmpz_mod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fmpz_mod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fmpz_randtest(c, state, 200);

            fmpz_mod_mpoly_add_fmpz(g, f, c, ctx);
            fmpz_mod_mpoly_assert_canonical(g, ctx);
            fmpz_mod_mpoly_sub_fmpz(h, g, c, ctx);
            fmpz_mod_mpoly_assert_canonical(h, ctx);
            if (!fmpz_mod_mpoly_equal(f, h, ctx))
            {
                flint_printf("FAIL: Check (f + c) - c = f\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }

            fmpz_mod_mpoly_set(g, f, ctx);
            fmpz_mod_mpoly_assert_canonical(f, ctx);
            fmpz_mod_mpoly_add_fmpz(f, f, c, ctx);
            fmpz_mod_mpoly_assert_canonical(f, ctx);
            fmpz_mod_mpoly_sub_fmpz(f, f, c, ctx);
            fmpz_mod_mpoly_assert_canonical(f, ctx);
            if (!fmpz_mod_mpoly_equal(f, g, ctx))
            {
                flint_printf("FAIL: Check aliasing\n");
                flint_printf("i = %wd, j = %wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_clear(c);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_clear(g, ctx);
        fmpz_mod_mpoly_clear(h, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
