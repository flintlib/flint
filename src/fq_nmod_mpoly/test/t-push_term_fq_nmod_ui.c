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

TEST_FUNCTION_START(fq_nmod_mpoly_push_term_fq_nmod_ui, state)
{
    slong i, j, k;

    /* Check pushback matches add */
    for (i = 0; i < 800 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f1, f2, m;
        flint_bitcnt_t exp_bits;
        ulong * exp, * exp2;
        slong len, nvars;
        fq_nmod_t c, c2;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f1, ctx);
        fq_nmod_mpoly_init(f2, ctx);
        fq_nmod_mpoly_init(m, ctx);
        fq_nmod_init(c, ctx->fqctx);
        fq_nmod_init(c2, ctx->fqctx);

        nvars = fq_nmod_mpoly_ctx_nvars(ctx);

        exp = (ulong *) flint_malloc(nvars*sizeof(ulong));
        exp2 = (ulong *) flint_malloc(nvars*sizeof(ulong));

        len = n_randint(state, 20);
        exp_bits = n_randint(state, FLINT_BITS) + 1;

        fq_nmod_mpoly_zero(f1, ctx);
        fq_nmod_mpoly_zero(f2, ctx);

        for (j = 0; j < len; j++)
        {
            /* get random term */
            fq_nmod_randtest(c, state, ctx->fqctx);
            for (k = 0; k < nvars; k++)
                exp[k] = n_randtest_bits(state, n_randint(state, exp_bits) + 1);

            /* add it to f1 */
            fq_nmod_mpoly_zero(m, ctx);
            fq_nmod_mpoly_set_coeff_fq_nmod_ui(m, c, exp, ctx);
            fq_nmod_mpoly_add(f1, f1, m, ctx);

            /* push it back on f2 */
            fq_nmod_mpoly_push_term_fq_nmod_ui(f2, c, exp, ctx);

            /* make sure last term matches */
            fq_nmod_mpoly_get_term_coeff_fq_nmod(c2, f2,
                                       fq_nmod_mpoly_length(f2, ctx) - 1, ctx);
            fq_nmod_mpoly_get_term_exp_ui(exp2, f2,
                                       fq_nmod_mpoly_length(f2, ctx) - 1, ctx);
            if (!fq_nmod_equal(c, c2, ctx->fqctx))
            {
                printf("FAIL\n");
                flint_printf("Check pushed coefficient matches\ni=%wd, j=%wd\n", i, j);
                fflush(stdout);
                flint_abort();
            }
            for (k = 0; k < nvars; k++)
            {
                if (exp[k] != exp2[k])
                {
                    printf("FAIL\n");
                    flint_printf("Check pushed exponent matches\ni=%wd, j=%wd\n", i, j);
                    fflush(stdout);
                    flint_abort();
                }
            }
        }

        fq_nmod_mpoly_sort_terms(f2, ctx);
        fq_nmod_mpoly_combine_like_terms(f2, ctx);
        fq_nmod_mpoly_assert_canonical(f2, ctx);

        if (!fq_nmod_mpoly_equal(f1, f2, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check pushback matches add\ni=%wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(c2, ctx->fqctx);
        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_mpoly_clear(f1, ctx);
        fq_nmod_mpoly_clear(f2, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);

        flint_free(exp2);
        flint_free(exp);
    }

    TEST_FUNCTION_END(state);
}
