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

TEST_FUNCTION_START(fmpz_mod_mpoly_get_set_str_pretty, state)
{
    slong i;

    {
        slong len1;
        flint_bitcnt_t exp_bits1;
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f, f1;
        char * str;
        const char * vars[] = {"x","y","z","w","u","v"};

        for (i = 0; i < flint_test_multiplier(); i++)
        {
            fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 6, 200);

            fmpz_mod_mpoly_init(f, ctx);
            fmpz_mod_mpoly_init(f1, ctx);

            for (len1 = 3; len1 < 1000; len1 += len1/2)
            {
                exp_bits1 = 200;
                fmpz_mod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);

                str = fmpz_mod_mpoly_get_str_pretty(f, vars, ctx);
                fmpz_mod_mpoly_set_str_pretty(f1, str, vars, ctx);
                flint_free(str);

                if (!fmpz_mod_mpoly_equal(f, f1, ctx))
                {
                    flint_printf("FAIL\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            fmpz_mod_mpoly_clear(f, ctx);
            fmpz_mod_mpoly_clear(f1, ctx);

            fmpz_mod_mpoly_ctx_clear(ctx);
        }
    }

    TEST_FUNCTION_END(state);
}
