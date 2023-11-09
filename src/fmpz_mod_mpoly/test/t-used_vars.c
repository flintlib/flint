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

TEST_FUNCTION_START(fmpz_mod_mpoly_used_vars, state)
{
    slong i, j, k;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fmpz_mod_mpoly_ctx_t ctx;
        fmpz_mod_mpoly_t f;
        fmpz_t one, fdeg;
        slong len, var;
        flint_bitcnt_t exp_bits;
        int * used;

        fmpz_mod_mpoly_ctx_init_rand_bits(ctx, state, 20, 200);
        fmpz_mod_mpoly_init(f, ctx);
        fmpz_init(fdeg);
        fmpz_init_set_ui(one, 1);
        used = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, int);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            len = n_randint(state, 200);
            exp_bits = n_randint(state, 150) + 2;

            fmpz_mod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

            for (k = n_randint(state, ctx->minfo->nvars); k > 0; k--)
            {
                var = n_randint(state, ctx->minfo->nvars);
                fmpz_mod_mpoly_evaluate_one_fmpz(f, f, var, one, ctx);
            }

            fmpz_mod_mpoly_used_vars(used, f, ctx);

            for (var = 0; var < ctx->minfo->nvars; var++)
            {
                fmpz_mod_mpoly_degree_fmpz(fdeg, f, var, ctx);
                if ((fmpz_sgn(fdeg) <= 0) != !used[var])
                {
                    flint_printf("FAIL: checked used matches degree\n");
                    flint_printf("var = %wd\n", var);
                    flint_printf("deg: "); fmpz_print(fdeg); flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        flint_free(used);
        fmpz_clear(one);
        fmpz_clear(fdeg);
        fmpz_mod_mpoly_clear(f, ctx);
        fmpz_mod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
