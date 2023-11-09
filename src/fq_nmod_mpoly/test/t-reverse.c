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

TEST_FUNCTION_START(fq_nmod_mpoly_reverse, state)
{
    int i, result;

    /* Check rev(rev(a)) == a */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        slong len, exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);

        len = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;

        fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        fq_nmod_mpoly_reverse(h, f, ctx);
        fq_nmod_mpoly_reverse(g, h, ctx);

        result = fq_nmod_mpoly_equal(f, g, ctx);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("Check rev(rev(a)) == a\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    /* Check aliasing */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g;
        slong len, exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);

        len = n_randint(state, 100);

        exp_bits = n_randint(state, 200) + 1;

        fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);

        fq_nmod_mpoly_reverse(g, f, ctx);
        fq_nmod_mpoly_reverse(g, g, ctx);

        result = fq_nmod_mpoly_equal(f, g, ctx);

        if (!result)
        {
            printf("FAIL\n");
            flint_printf("Check aliasing\ni = %wd\n", i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
