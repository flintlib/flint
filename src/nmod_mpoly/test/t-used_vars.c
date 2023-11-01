/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mpoly.h"

TEST_FUNCTION_START(nmod_mpoly_used_vars, state)
{
    slong i, j, k;

    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        nmod_mpoly_ctx_t ctx;
        nmod_mpoly_t f;
        fmpz_t fdeg;
        slong len, var;
        flint_bitcnt_t exp_bits;
        int * used;

        nmod_mpoly_ctx_init_rand(ctx, state, 10, n_randint(state, 100) + 1);
        nmod_mpoly_init(f, ctx);
        fmpz_init(fdeg);
        used = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, int);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            len = n_randint(state, 500);
            exp_bits = n_randint(state, 20) + 2;

            nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

            for (k = n_randint(state, ctx->minfo->nvars); k > 0; k--)
            {
                var = n_randint(state, ctx->minfo->nvars);
                nmod_mpoly_evaluate_one_ui(f, f, var, 1, ctx);
            }

            nmod_mpoly_used_vars(used, f, ctx);

            for (var = 0; var < ctx->minfo->nvars; var++)
            {
                nmod_mpoly_degree_fmpz(fdeg, f, var, ctx);
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
        fmpz_clear(fdeg);
        nmod_mpoly_clear(f, ctx);
        nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
