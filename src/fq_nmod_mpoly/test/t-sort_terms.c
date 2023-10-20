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

TEST_FUNCTION_START(fq_nmod_mpoly_sort_terms, state)
{
    int i, j, result;

    /* Check scramble and sort */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f, g;
        slong len;
        flint_bitcnt_t exp_bits;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 20);

        fq_nmod_mpoly_init(f, ctx);
        fq_nmod_mpoly_init(g, ctx);

        len = n_randint(state, 200);
        exp_bits = n_randint(state, 200) + 1;

        for (j = 0; j < 4; j++)
        {
            slong d = fq_nmod_ctx_degree(ctx->fqctx);
            slong N, k;

            fq_nmod_mpoly_randtest_bits(f, state, len, exp_bits, ctx);
            fq_nmod_mpoly_set(g, f, ctx);

            N = mpoly_words_per_exp(f->bits, ctx->minfo);
            for (k = WORD(0); k < f->length; k++)
            {
                ulong a, b;
                a = n_randint(state, f->length);
                b = n_randint(state, f->length);
                _n_fq_swap(f->coeffs + d*a, f->coeffs + d*b, d);
                mpoly_monomial_swap(f->exps + N*a, f->exps + N*b, N);
            }

            fq_nmod_mpoly_sort_terms(f, ctx);
            result = fq_nmod_mpoly_equal(f, g, ctx);
            if (!result)
            {
                printf("FAIL\n");
                flint_printf("Check scramble and sort\ni = %wd, j = %wd\n", i ,j);
                fflush(stdout);
                flint_abort();
            }
        }

        fq_nmod_mpoly_clear(f, ctx);
        fq_nmod_mpoly_clear(g, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
