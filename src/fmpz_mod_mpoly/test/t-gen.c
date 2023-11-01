/*
    Copyright (C) 2022 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_gen, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f1, f2;
        slong nvars, len, exp_bits, k1, k2;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        nvars = fmpz_mod_mpoly_ctx_nvars(ctx);
        if (nvars < 1)
        {
            fmpz_mod_mpoly_ctx_clear(ctx);
            continue;
        }

        fmpz_mod_mpoly_init(f1, ctx);
        fmpz_mod_mpoly_init(f2, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;

        fmpz_mod_mpoly_randtest_bits(f1, state, len, exp_bits, ctx);
        fmpz_mod_mpoly_randtest_bits(f2, state, len, exp_bits, ctx);

        k1 = n_randint(state, nvars);
        k2 = n_randint(state, nvars);
        fmpz_mod_mpoly_gen(f1, k1, ctx);
        fmpz_mod_mpoly_assert_canonical(f1, ctx);
        fmpz_mod_mpoly_gen(f2, k2, ctx);
        fmpz_mod_mpoly_assert_canonical(f2, ctx);

        result = 1;
        result = result && fmpz_mod_mpoly_is_gen(f1, k1, ctx);
        result = result && fmpz_mod_mpoly_is_gen(f1, -WORD(1), ctx);
        result = result && fmpz_mod_mpoly_is_gen(f2, k2, ctx);
        result = result && fmpz_mod_mpoly_is_gen(f2, -WORD(1), ctx);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("Check one generator\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_mul(f1, f1, f2, ctx);
        result = 0;
        result += !!fmpz_mod_mpoly_is_gen(f1, k1, ctx);
        result += !!fmpz_mod_mpoly_is_gen(f1, k2, ctx);
        result += !!fmpz_mod_mpoly_is_gen(f1, -WORD(1), ctx);

        if (result != (fmpz_is_one(fmpz_mod_mpoly_ctx_modulus(ctx)) ? 3 : 0))
        {
            printf("FAIL\n");
            flint_printf("Check product of two generators\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fmpz_mod_mpoly_clear(f1, ctx);
        fmpz_mod_mpoly_clear(f2, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
