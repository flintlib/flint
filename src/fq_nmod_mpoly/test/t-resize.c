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

TEST_FUNCTION_START(fq_nmod_mpoly_resize, state)
{
    slong i, j, k;

    /* Check pushback matches add */
    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fq_nmod_mpoly_ctx_t ctx;
        fq_nmod_mpoly_t f1, f2, f3, m;
        flint_bitcnt_t exp_bits;
        ulong * exp, * exp2;
        slong len, nvars;
        fq_nmod_t c;

        fq_nmod_mpoly_ctx_init_rand(ctx, state, 10, FLINT_BITS, 10);
        fq_nmod_mpoly_init(f1, ctx);
        fq_nmod_mpoly_init(f2, ctx);
        fq_nmod_mpoly_init(f3, ctx);
        fq_nmod_mpoly_init(m, ctx);
        fq_nmod_init(c, ctx->fqctx);

        nvars = fq_nmod_mpoly_ctx_nvars(ctx);

        exp = (ulong *) flint_malloc(nvars*sizeof(ulong));
        exp2 = (ulong *) flint_malloc(nvars*sizeof(ulong));

        len = n_randint(state, 200);
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

            /* manually push it on f3 */
            fq_nmod_mpoly_resize(f3, j + 1 + n_randint(state, 10), ctx);
            fq_nmod_mpoly_set_term_coeff_fq_nmod(f3, j, c, ctx);
            fq_nmod_mpoly_set_term_exp_ui(f3, j, exp, ctx);
        }

        fq_nmod_mpoly_sort_terms(f2, ctx);
        fq_nmod_mpoly_combine_like_terms(f2, ctx);
        fq_nmod_mpoly_assert_canonical(f2, ctx);

        fq_nmod_mpoly_sort_terms(f3, ctx);
        fq_nmod_mpoly_combine_like_terms(f3, ctx);
        fq_nmod_mpoly_assert_canonical(f3, ctx);

        if (!fq_nmod_mpoly_equal(f1, f2, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check pushback matches add\ni=%wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        if (!fq_nmod_mpoly_equal(f1, f3, ctx))
        {
            printf("FAIL\n");
            flint_printf("Check resize+setterm matches add\ni=%wd\n",i);
            fflush(stdout);
            flint_abort();
        }

        fq_nmod_clear(c, ctx->fqctx);
        fq_nmod_mpoly_clear(f1, ctx);
        fq_nmod_mpoly_clear(f2, ctx);
        fq_nmod_mpoly_clear(f3, ctx);
        fq_nmod_mpoly_clear(m, ctx);
        fq_nmod_mpoly_ctx_clear(ctx);

        flint_free(exp2);
        flint_free(exp);
    }

    TEST_FUNCTION_END(state);
}
