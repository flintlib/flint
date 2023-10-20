/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_used_vars, state)
{
    slong i, j, k;

    for (i = 0; i < 20 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f;
        fmpz_t fdeg;
        slong len, var;
        fq_nmod_t one;
        flint_bitcnt_t exp_bits;
        int * used;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f, ctx);
        fmpz_init(fdeg);
        fq_nmod_init(one, ctx->fqctx);
        fq_nmod_sub_one(one, one, ctx->fqctx);
        fq_nmod_neg(one, one, ctx->fqctx);
        used = FLINT_ARRAY_ALLOC(ctx->minfo->nvars, int);

        for (j = 0; j < ctx->minfo->nvars; j++)
        {
            len = n_randint(state, 200);
            exp_bits = n_randint(state, 200) + 2;

            fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

            for (k = n_randint(state, ctx->minfo->nvars); k > 0; k--)
            {
                var = n_randint(state, ctx->minfo->nvars);
                fq_nmod_mpoly_evaluate_one_fq_nmod(f, f, var, one, ctx);
            }

            fq_nmod_mpoly_used_vars(used, f, ctx);

            for (var = 0; var < ctx->minfo->nvars; var++)
            {
                fq_nmod_mpoly_degree_fmpz(fdeg, f, var, ctx);
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
        fq_nmod_clear(one, ctx->fqctx);
        fmpz_clear(fdeg);
        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
