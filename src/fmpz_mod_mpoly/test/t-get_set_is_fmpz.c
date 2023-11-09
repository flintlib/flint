/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz_mod_mpoly.h"

TEST_FUNCTION_START(fmpz_mod_mpoly_get_set_is_fmpz, state)
{
    slong i;

    /* Set to random integer and compare */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f;
        fmpz_t c, d;

        fmpz_init(c);
        fmpz_init(d);

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        fmpz_mod_mpoly_init(f, ctx);

        fmpz_mod_mpoly_randtest_bits(f, state, n_randint(state, 10),
                                                   n_randint(state, 200), ctx);

        if (fmpz_mod_mpoly_is_fmpz(f, ctx))
        {
            fmpz_mod_mpoly_get_fmpz(c, f, ctx);
            fmpz_add(c, c, fmpz_mod_mpoly_ctx_modulus(ctx));
            if (!fmpz_mod_mpoly_equal_fmpz(f, c, ctx))
            {
                flint_printf("FAIL: Check is_fmpz and get_fmpz match\n");
                fflush(stdout);
                flint_abort();
            }
        }

        fmpz_randtest(c, state, n_randint(state, 300));
        fmpz_mod_mpoly_set_fmpz(f, c, ctx);
        if (!fmpz_mod_mpoly_is_fmpz(f, ctx))
        {
            flint_printf("FAIL: Check set_fmpz makes is_fmpz true\n");
            fflush(stdout);
            flint_abort();
        }
        fmpz_mod_mpoly_get_fmpz(d, f, ctx);
        if (!fmpz_mod_equal_fmpz(d, c, ctx->ffinfo))
        {
            flint_printf("FAIL: Check get_fmpz matches set_fmpz\n");
            fflush(stdout);
            flint_abort();
        }

        fmpz_clear(c);
        fmpz_clear(d);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
