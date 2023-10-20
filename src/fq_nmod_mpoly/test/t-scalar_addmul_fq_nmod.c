/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_scalar_addmul_fq_nmod, state)
{
    slong i, j;

    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h, k;
        fq_nmod_t c;
        slong len, len1, len2;
        slong exp_bits, exp_bits1, exp_bits2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_mpoly_init(k, ctx);

        fq_nmod_init(c, ctx->fqctx);

        len = n_randint(state, 100);
        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 2;
        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;

        for (j = 0; j < 10; j++)
        {
            fq_nmod_randtest(c, state, ctx->fqctx);
            fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(k, state, len, exp_bits, ctx);

            fq_nmod_mpoly_scalar_addmul_fq_nmod(f, g, h, c, ctx);

            fq_nmod_mpoly_scalar_mul_fq_nmod(k, h, c, ctx);
            fq_nmod_mpoly_add(k, k, g, ctx);

            if (!fq_nmod_mpoly_equal(f, k, ctx))
            {
                flint_printf("FAIL: check definition\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_set(k, g, ctx);
            fq_nmod_mpoly_scalar_addmul_fq_nmod(k, k, h, c, ctx);
            if (!fq_nmod_mpoly_equal(f, k, ctx))
            {
                flint_printf("FAIL: check aliasing first argument\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }

            fq_nmod_mpoly_set(k, h, ctx);
            fq_nmod_mpoly_scalar_addmul_fq_nmod(k, g, k, c, ctx);
            if (!fq_nmod_mpoly_equal(f, k, ctx))
            {
                flint_printf("FAIL: check aliasing second argument\n");
                flint_printf("i = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_clear(c, ctx->fqctx);

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_clear(k, ctx);

        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
