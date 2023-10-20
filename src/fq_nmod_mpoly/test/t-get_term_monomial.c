/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fq_nmod_mpoly.h"

TEST_FUNCTION_START(fq_nmod_mpoly_get_term_monomial, state)
{
    int i, j;

    /* Check getting a coeff by its monomial */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_t c, d;
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g, h;
        flint_bitcnt_t exp_bits1, exp_bits2, exp_bits3;
        slong len1, len2, len3;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 20, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);
        fq_nmod_mpoly_init(h, ctx);
        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_init(d, ctx->fqctx);

        len1 = n_randint(state, 100);
        len2 = n_randint(state, 100);
        len3 = n_randint(state, 100);

        exp_bits1 = n_randint(state, 200) + 2;
        exp_bits2 = n_randint(state, 200) + 2;
        exp_bits3 = n_randint(state, 200) + 2;

        fq_nmod_mpoly_randtest_bits(f, state, len1, exp_bits1, ctx);
        fq_nmod_mpoly_randtest_bits(g, state, len2, exp_bits2, ctx);
        fq_nmod_mpoly_randtest_bits(h, state, len3, exp_bits3, ctx);

        fq_nmod_mpoly_repack_bits(h, f,
                                f->bits + n_randint(state, 2*FLINT_BITS), ctx);

        for (j = fq_nmod_mpoly_length(f, ctx) - 1; j >= 0; j--)
        {
            fq_nmod_mpoly_get_term_monomial(g, f, j, ctx);
            fq_nmod_mpoly_repack_bits(g, g,
                                  g->bits + n_randint(state, FLINT_BITS), ctx);
            fq_nmod_mpoly_get_term_coeff_fq_nmod(d, f, j, ctx);
            fq_nmod_mpoly_get_coeff_fq_nmod_monomial(c, h, g, ctx);

            if (!fq_nmod_equal(c, d, ctx->fqctx))
            {
                flint_printf("FAIL\nCheck getting a coeff by its monomial\ni = %wd\n", i);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_clear(h, ctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_clear(d, ctx->fqctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
