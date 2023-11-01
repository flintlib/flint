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

TEST_FUNCTION_START(fq_nmod_mpoly_add_sub_fq_nmod, state)
{
    int i, j, result;

    /* Check (f + a) - a = f */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        fq_nmod_t c;
        slong len, exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_init(c, ctx->fqctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;

        for (j = 0; j < 10; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(g, state, len, exp_bits, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            fq_nmod_randtest(c, state, ctx->fqctx);

            fq_nmod_mpoly_add_fq_nmod(g, f, c, ctx);
            fq_nmod_mpoly_assert_canonical(g, ctx);

            fq_nmod_mpoly_sub_fq_nmod(h, g, c, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            result = fq_nmod_mpoly_equal(f, h, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check (f + a) - a = f\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 10 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        fq_nmod_t c;
        slong len, exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_init(c, ctx->fqctx);

        len = n_randint(state, 100);
        exp_bits = n_randint(state, 200) + 1;

        for (j = 0; j < 10; j++)
        {
            fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fq_nmod_mpoly_set(g, f, ctx);
            fq_nmod_mpoly_randtest_bits(h, state, len, exp_bits, ctx);

            fq_nmod_randtest(c, state, ctx->fqctx);

            fq_nmod_mpoly_set_fq_nmod(h, c, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_add(h, f, h, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            fq_nmod_mpoly_add_fq_nmod(f, f, c, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);

            result = fq_nmod_mpoly_equal(f, h, ctx);

            fq_nmod_mpoly_set_fq_nmod(h, c, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);
            fq_nmod_mpoly_sub(h, f, h, ctx);
            fq_nmod_mpoly_assert_canonical(h, ctx);

            fq_nmod_mpoly_sub_fq_nmod(f, f, c, ctx);
            fq_nmod_mpoly_assert_canonical(f, ctx);

            result = result && fq_nmod_mpoly_equal(f, h, ctx);

            result = result && fq_nmod_mpoly_equal(f, g, ctx);

            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check aliasing\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
