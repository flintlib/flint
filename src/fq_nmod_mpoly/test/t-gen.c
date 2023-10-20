/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_gen, state)
{
    int i, result;

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f1, f2;
        slong len, k1, k2;
        flint_bitcnt_t exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        if (fq_nmod_mpoly_ctx_nvars(ctx) < 1)
        {
            fq_nmod_mpoly_ctx_clear(ctx);
            continue;
        }

        fq_nmod_mpoly_init(f1, ctx);
        fq_nmod_mpoly_init(f2, ctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;

        fq_nmod_mpoly_randtest_bits(f1, state, len, exp_bits, ctx);
        fq_nmod_mpoly_randtest_bits(f2, state, len, exp_bits, ctx);

        k1 = n_randint(state, ctx->minfo->nvars);
        k2 = n_randint(state, ctx->minfo->nvars);
        fq_nmod_mpoly_gen(f1, k1, ctx);
        fq_nmod_mpoly_assert_canonical(f1, ctx);
        fq_nmod_mpoly_gen(f2, k2, ctx);
        fq_nmod_mpoly_assert_canonical(f2, ctx);

        result = 1;
        result = result && fq_nmod_mpoly_is_gen(f1, k1, ctx);
        result = result && fq_nmod_mpoly_is_gen(f1, -WORD(1), ctx);
        result = result && fq_nmod_mpoly_is_gen(f2, k2, ctx);
        result = result && fq_nmod_mpoly_is_gen(f2, -WORD(1), ctx);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("Check one generator\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_mul(f1, f1, f2, ctx);
        result = 1;
        result = result && !fq_nmod_mpoly_is_gen(f1, k1, ctx);
        result = result && !fq_nmod_mpoly_is_gen(f1, k2, ctx);
        result = result && !fq_nmod_mpoly_is_gen(f1, -WORD(1), ctx);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("Check product of two generators\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_clear(f1, ctx);
        fq_nmod_mpoly_clear(f2, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
